/*
 * Multi-threaded mpileup, enabled with --threads.
 *
 * The regions are split into chunks which are processed independently by
 * worker threads, each with its own file handles, pileup iterator and calling
 * structures, and the main thread outputs the chunks in order.
 *
 * The output is byte-identical to a single-threaded run. A worker starts the
 * pileup at a point before its chunk (proc_beg), processes the columns before
 * the chunk without outputting them, and so every read which reaches the chunk
 * goes through exactly the same history as in a sequential run (BAQ, quality
 * adjustment of overlapping mates, -d subsampling). The start point is chosen
 * so that
 *   - reads which start before proc_beg, and their overlapping mates, end
 *     before the chunk begins
 *   - there are fewer reads around proc_beg than the -d limit, so that no read
 *     can be dropped there by the subsampling and the set of reads in the
 *     pileup is the same as in a sequential run.
 * At the other end, the pileup reads ahead one read past the current column,
 * which can change the qualities (overlapping mates) and therefore BAQ of the
 * reads in the column. Therefore reads are fetched beyond the end of the chunk
 * by more than the longest read, except at the end of a region where the
 * sequential run stops as well.
 * If no such point can be found (e.g. the depth is constantly above -d), the
 * chunk is merged with the previous one.
 *
 * Genotype likelihoods at sites with more than 255 bases are calculated from a
 * random subset of the bases (errmod_cal), so the output depends on how many
 * random numbers were drawn before. The number of draws in a chunk does not
 * depend on the random numbers, so each chunk counts its draws and, when all
 * the preceding chunks are done, the state of the generator at the start of
 * the chunk is known (hts_rand48_skip). A chunk which drew random numbers
 * before that state was known is then processed again.
 */

typedef struct
{
    int tid;
    hts_pos_t reg_beg;      // start of the enclosing region, processing never starts before it
    hts_pos_t proc_beg;     // where the pileup starts
    hts_pos_t beg, end;     // records are output for columns in [beg,end]
    hts_pos_t span;         // the longest read seen around the chunk start
    int is_last;            // the last chunk of a region: reads are fetched up to the end only
}
mplp_chunk_t;

typedef struct
{
    mplp_conf_t *conf;
    bam_hdr_t *hdr;
    mplp_chunk_t *chunk;
    int nchunk, mchunk;
    mplp_chunk_res_t *res;  // one per chunk
    int *status;            // per chunk: 0 pending, 1 done, -1 failed
    int *final;             // per chunk: the output does not need to be regenerated
    uint64_t *ndraw;        // per chunk: the number of random numbers drawn
    unsigned short (*rng)[3];   // per chunk: generator state at the chunk start, valid for chunks < nrng
    int nrng;
    int *rerun, nrerun;     // chunks to process again with the correct generator state
    int finished;
    kstring_t *pool;        // output buffers for reuse
    int npool, mpool;
    int next, nwritten, max_inflight, failed;
    pthread_mutex_t lock;
    pthread_cond_t cond;
}
mplp_par_t;

typedef struct
{
    mplp_par_t *par;
    mplp_conf_t conf;
    mplp_ref_t ref;
    pthread_t thread;
}
mplp_worker_t;

#define MPLP_PAR_CHUNK_SIZE 100000

static int mpileup_parallel_ok(mplp_conf_t *conf, int nregs)
{
    int i;
    if ( conf->n_threads < 1 ) return 0;
    if ( conf->gvcf || conf->write_index ) return 0;
    if ( conf->reg && !nregs ) return 0;
    unsigned short xseed[3];
    if ( hts_drand48_get_state(xseed) < 0 ) return 0;

    htsFile *fp = conf->bcf_fp;
    if ( fp->format.format==bcf ) conf->par_is_vcf = 0;
    else if ( fp->format.compression==no_compression ) conf->par_is_vcf = 1;
    else if ( fp->format.compression==bgzf ) conf->par_is_vcf = 2;
    else return 0;

    // all input files must be indexed
    for (i=0; i<conf->nfiles; i++)
    {
        hts_idx_t *idx = sam_index_load(conf->mplp_data[i]->fp, conf->files[i]);
        if ( !idx )
        {
            fprintf(stderr,"[%s] Note: %s is not indexed, running single-threaded\n", __func__, conf->files[i]);
            return 0;
        }
        hts_idx_destroy(idx);
    }
    return 1;
}

// Count reads overlapping [beg,end] and report the longest reference span seen
static int mplp_par_count(samFile *fp, hts_idx_t *idx, bam1_t *b, int tid, hts_pos_t beg, hts_pos_t end, hts_pos_t *max_span)
{
    hts_itr_t *itr = sam_itr_queryi(idx, tid, beg < 0 ? 0 : beg, end+1);
    if ( !itr ) error("[%s] Error: failed to query the index\n", __func__);
    int n = 0, ret;
    while ( (ret=sam_itr_next(fp, itr, b))>=0 )
    {
        if ( b->core.flag & BAM_FUNMAP ) continue;
        hts_pos_t span = bam_endpos(b) - b->core.pos;
        if ( *max_span < span ) *max_span = span;
        n++;
    }
    hts_itr_destroy(itr);
    if ( ret < -1 ) error("[%s] Error: failed to read the input\n", __func__);
    return n;
}

// Find where to start the pileup for a chunk beginning at beg, returns -1 if
// no suitable point was found
static int mplp_par_find_start(mplp_conf_t *conf, samFile **fp, hts_idx_t **idx, bam1_t *b,
        int tid, hts_pos_t reg_beg, hts_pos_t beg, hts_pos_t *proc_beg, hts_pos_t *max_read_span)
{
    hts_pos_t span = 0;     // the longest read around the start point
    int itry = 0, i;
    while ( itry < 16 )
    {
        hts_pos_t p = beg - 2*span - 100 - itry*(span+100);
        if ( p <= reg_beg ) { *proc_beg = reg_beg; return 0; }
        int quiet = 1;
        hts_pos_t max_span = span;
        for (i=0; i<conf->nfiles; i++)
        {
            int n = mplp_par_count(fp[i], idx[i], b, tid, p - span, p + span, &max_span);
            if ( (int64_t)n + 2 > conf->max_depth ) quiet = 0;
        }
        if ( max_span > span ) { span = max_span; continue; }   // re-check with the longer span
        *max_read_span = span;
        if ( quiet ) { *proc_beg = p; return 0; }
        itry++;
    }
    return -1;
}

static void mplp_par_add_chunks(mplp_par_t *par, samFile **fp, hts_idx_t **idx, bam1_t *b,
        int tid, hts_pos_t reg_beg, hts_pos_t reg_end, hts_pos_t emit_end)
{
    mplp_conf_t *conf = par->conf;
    hts_pos_t beg = reg_beg;
    while ( 1 )
    {
        hts_pos_t end = beg + conf->par_chunk_size - 1;
        if ( end >= reg_end ) end = emit_end;

        hts_pos_t proc_beg = beg, span = 0;
        if ( beg > reg_beg && mplp_par_find_start(conf, fp, idx, b, tid, reg_beg, beg, &proc_beg, &span)<0 )
        {
            // no safe starting point, extend the previous chunk
            par->chunk[par->nchunk-1].end = end;
            par->chunk[par->nchunk-1].is_last = end==emit_end ? 1 : 0;
        }
        else
        {
            hts_expand(mplp_chunk_t, par->nchunk+1, par->mchunk, par->chunk);
            mplp_chunk_t *chunk = &par->chunk[par->nchunk++];
            chunk->tid = tid;
            chunk->reg_beg  = reg_beg;
            chunk->proc_beg = proc_beg;
            chunk->beg = beg;
            chunk->end = end;
            chunk->span = span;
            chunk->is_last = end==emit_end ? 1 : 0;
        }
        if ( end==emit_end ) break;
        beg = end + 1;
    }
}

static void mplp_par_plan(mplp_par_t *par, int nregs)
{
    mplp_conf_t *conf = par->conf;
    samFile **fp = (samFile**) calloc(conf->nfiles, sizeof(*fp));
    hts_idx_t **idx = (hts_idx_t**) calloc(conf->nfiles, sizeof(*idx));
    bam1_t *b = bam_init1();
    int i;
    for (i=0; i<conf->nfiles; i++)
    {
        fp[i] = sam_open(conf->files[i], "rb");
        if ( !fp[i] ) error("[%s] Error: failed to open %s\n", __func__, conf->files[i]);
        if ( conf->fai_fname ) hts_set_fai_filename(fp[i], conf->fai_fname);
        bam_hdr_t *h = sam_hdr_read(fp[i]);
        if ( !h ) error("[%s] Error: failed to read the header of %s\n", __func__, conf->files[i]);
        bam_hdr_destroy(h);
        idx[i] = sam_index_load(fp[i], conf->files[i]);
        if ( !idx[i] ) error("[%s] Error: failed to load the index of %s\n", __func__, conf->files[i]);
    }

    if ( nregs )
    {
        // the region iterator is already positioned at the first region
        do
        {
            int tid = sam_hdr_name2tid(par->hdr, conf->reg_itr->seq);
            if ( tid < 0 )
                error("[E::%s] the sequence \"%s\" not found: %s\n",__func__,conf->reg_itr->seq,conf->files[0]);
            mplp_par_add_chunks(par, fp, idx, b, tid, conf->reg_itr->beg, conf->reg_itr->end, conf->reg_itr->end);
        }
        while ( regitr_loop(conf->reg_itr) );
    }
    else
    {
        // positions beyond the end of the sequence are output in the sequential mode too
        for (i=0; i<par->hdr->n_targets; i++)
        {
            // skip sequences with no reads
            int j, has_reads = 0;
            for (j=0; j<conf->nfiles && !has_reads; j++)
            {
                hts_itr_t *itr = sam_itr_queryi(idx[j], i, 0, HTS_POS_MAX);
                if ( !itr ) error("[%s] Error: failed to query the index\n", __func__);
                if ( sam_itr_next(fp[j], itr, b) >= 0 ) has_reads = 1;
                hts_itr_destroy(itr);
            }
            if ( !has_reads ) continue;
            mplp_par_add_chunks(par, fp, idx, b, i, 0, par->hdr->target_len[i] ? par->hdr->target_len[i]-1 : 0, HTS_POS_MAX);
        }
    }

    bam_destroy1(b);
    for (i=0; i<conf->nfiles; i++)
    {
        hts_idx_destroy(idx[i]);
        sam_close(fp[i]);
    }
    free(idx);
    free(fp);
}

static void mplp_par_worker_init(mplp_worker_t *w, mplp_par_t *par)
{
    mplp_conf_t *main_conf = par->conf, *conf = &w->conf;
    *conf = *main_conf;
    w->par = par;
    mplp_ref_t ref = MPLP_REF_INIT;
    w->ref = ref;

    conf->progress = NULL;
    conf->iter = NULL;
    conf->reg = NULL;
    conf->reg_itr = NULL;
    conf->par_res = NULL;
    memset(&conf->buf, 0, sizeof(conf->buf));
    memset(&conf->bc, 0, sizeof(conf->bc));
    conf->bed_itr = conf->bed ? regitr_init(conf->bed) : NULL;
    // a private copy, htslib caches data in the header when formatting VCF
    conf->bcf_hdr = bcf_hdr_dup(main_conf->bcf_hdr);
    if ( !conf->bcf_hdr ) error("[%s] Error: failed to duplicate the header\n", __func__);
    conf->fai = NULL;
    if ( conf->fai_fname && !(conf->fai = fai_load(conf->fai_fname)) )
        error("[%s] Error: failed to load the reference %s\n", __func__, conf->fai_fname);

    int i;
    conf->mplp_data = (mplp_aux_t**) calloc(conf->nfiles, sizeof(mplp_aux_t*));
    conf->plp = (const bam_pileup1_t**) calloc(conf->nfiles, sizeof(bam_pileup1_t*));
    conf->n_plp = (int*) calloc(conf->nfiles, sizeof(int));
    for (i=0; i<conf->nfiles; i++)
    {
        mplp_aux_t *ma = conf->mplp_data[i] = (mplp_aux_t*) calloc(1, sizeof(mplp_aux_t));
        ma->fp = sam_open(conf->files[i], "rb");
        if ( !ma->fp ) error("[%s] Error: failed to open %s\n", __func__, conf->files[i]);
        if ( hts_set_opt(ma->fp, CRAM_OPT_DECODE_MD, 1) ) error("Failed to set CRAM_OPT_DECODE_MD value\n");
        if ( conf->fai_fname && hts_set_fai_filename(ma->fp, conf->fai_fname)!=0 )
            error("[%s] Error: failed to process %s\n", __func__, conf->fai_fname);
        bam_hdr_t *h = sam_hdr_read(ma->fp);
        if ( !h ) error("[%s] Error: failed to read the header of %s\n", __func__, conf->files[i]);
        bam_hdr_destroy(h);
        ma->h = par->hdr;
        ma->conf = conf;
        ma->ref = &w->ref;
        ma->bam_id = main_conf->mplp_data[i]->bam_id;
        ma->idx = sam_index_load(ma->fp, conf->files[i]);
        if ( !ma->idx ) error("[%s] Error: failed to load the index of %s\n", __func__, conf->files[i]);
    }

    conf->gplp = (mplp_pileup_t *) calloc(1,sizeof(mplp_pileup_t));
    conf->gplp->n = main_conf->gplp->n;
    conf->gplp->n_plp = (int*) calloc(conf->gplp->n, sizeof(int));
    conf->gplp->m_plp = (int*) calloc(conf->gplp->n, sizeof(int));
    conf->gplp->plp = (bam_pileup1_t**) calloc(conf->gplp->n, sizeof(bam_pileup1_t*));

    mplp_init_call(conf, main_conf->bc.n);
    conf->bcf_rec = bcf_init1();
    conf->par_rec = bcf_init1();

    // The error model tables are large (32MB) and read-only, share them
    errmod_destroy(conf->bca->e);
    conf->bca->e = main_conf->bca->e;
}

static void mplp_par_worker_destroy(mplp_worker_t *w)
{
    mplp_conf_t *conf = &w->conf;
    int i;
    bcf_destroy1(conf->bcf_rec);
    bcf_destroy1(conf->par_rec);
    conf->bca->e = NULL;    // shared with the main thread
    mplp_destroy_call(conf);
    for (i = 0; i < conf->gplp->n; ++i) free(conf->gplp->plp[i]);
    free(conf->gplp->plp); free(conf->gplp->n_plp); free(conf->gplp->m_plp); free(conf->gplp);
    for (i=0; i<conf->nfiles; i++)
    {
        mplp_aux_t *ma = conf->mplp_data[i];
        if ( ma->iter ) hts_itr_destroy(ma->iter);
        hts_idx_destroy(ma->idx);
        sam_close(ma->fp);
        free(ma);
    }
    free(conf->mplp_data); free(conf->plp); free(conf->n_plp);
    free(conf->buf.s);
    if ( conf->bed_itr ) regitr_destroy(conf->bed_itr);
    if ( conf->fai ) fai_destroy(conf->fai);
    bcf_hdr_destroy(conf->bcf_hdr);
    free(w->ref.ref[0]);
    free(w->ref.ref[1]);
}

static int mplp_par_run_chunk(mplp_worker_t *w, mplp_chunk_t *chunk)
{
    mplp_conf_t *conf = &w->conf;
    mplp_chunk_res_t *res = conf->par_res;
    hts_pos_t proc_beg = chunk->proc_beg;
    hts_pos_t lookahead = 2*chunk->span > 1000 ? 2*chunk->span : 1000;
    int i;
    while ( 1 )
    {
        res->buf.l = 0;
        res->nlen = 0;
        hts_drand48_set_thread_state(conf->par_rng_cur, NULL);  // until the first output column
        hts_pos_t itr_end;
        if ( chunk->end==HTS_POS_MAX ) itr_end = HTS_POS_MAX;
        else if ( chunk->is_last ) itr_end = chunk->end + 1;
        else itr_end = chunk->end + 1 + lookahead;
        for (i=0; i<conf->nfiles; i++)
        {
            mplp_aux_t *ma = conf->mplp_data[i];
            if ( ma->iter ) hts_itr_destroy(ma->iter);
            ma->iter = sam_itr_queryi(ma->idx, chunk->tid, proc_beg, itr_end);
            if ( !ma->iter ) return -1;
            ma->track = 1;
            ma->proc_beg = proc_beg;
            ma->taint_end = -1;
            ma->max_span = 0;
        }
        conf->iter = bam_mplp_init(conf->nfiles, mplp_func, (void**)conf->mplp_data);
        if ( conf->flag & MPLP_SMART_OVERLAPS ) bam_mplp_init_overlaps(conf->iter);
        bam_mplp_set_maxcnt(conf->iter, conf->max_depth);
        bam_mplp_constructor(conf->iter, pileup_constructor);
        bam_mplp_destructor(conf->iter, pileup_destructor);

        int ret = mpileup_reg(conf, proc_beg, chunk->end, chunk->beg);

        bam_mplp_destroy(conf->iter);
        conf->iter = NULL;
        if ( ret < 0 ) return -1;

        // Verify that the reads outside of the processed range could not affect the output
        hts_pos_t taint_end = -1, max_span = 0;
        for (i=0; i<conf->nfiles; i++)
        {
            if ( taint_end < conf->mplp_data[i]->taint_end ) taint_end = conf->mplp_data[i]->taint_end;
            if ( max_span < conf->mplp_data[i]->max_span ) max_span = conf->mplp_data[i]->max_span;
        }
        int ok = 1;
        if ( max_span >= lookahead && itr_end!=HTS_POS_MAX && !chunk->is_last )
        {
            lookahead = 2*max_span;
            ok = 0;
        }
        if ( proc_beg > chunk->reg_beg && taint_end >= 0 && taint_end + max_span > chunk->beg )
        {
            fprintf(stderr,"[%s] Note: re-running %s:%"PRIhts_pos" from the region start\n", __func__,
                    w->par->hdr->target_name[chunk->tid], chunk->beg+1);
            proc_beg = chunk->reg_beg;
            ok = 0;
        }
        if ( ok ) return 0;
    }
}

static void *mplp_par_worker(void *arg)
{
    mplp_worker_t *w = (mplp_worker_t*) arg;
    mplp_par_t *par = w->par;

    // never use the shared global state of the random number generator
    memset(w->conf.par_rng_cur, 0, sizeof(w->conf.par_rng_cur));
    hts_drand48_set_thread_state(w->conf.par_rng_cur, NULL);

    while ( 1 )
    {
        pthread_mutex_lock(&par->lock);
        while ( !par->failed && !par->finished && !par->nrerun
                && !(par->next < par->nchunk && par->next < par->nwritten + par->max_inflight) )
            pthread_cond_wait(&par->cond, &par->lock);
        if ( par->failed || (par->finished && !par->nrerun) )
        {
            pthread_mutex_unlock(&par->lock);
            break;
        }
        int ichunk;
        if ( par->nrerun ) ichunk = par->rerun[--par->nrerun];
        else
        {
            ichunk = par->next++;
            if ( par->npool ) par->res[ichunk].buf = par->pool[--par->npool];
        }
        int exact = ichunk < par->nrng ? 1 : 0;
        if ( exact ) memcpy(w->conf.par_rng, par->rng[ichunk], sizeof(w->conf.par_rng));
        else memset(w->conf.par_rng, 0, sizeof(w->conf.par_rng));
        pthread_mutex_unlock(&par->lock);

        w->conf.par_res = &par->res[ichunk];
        w->conf.par_ndraw = 0;
        int ret = mplp_par_run_chunk(w, &par->chunk[ichunk]);
        w->conf.par_res = NULL;
        hts_drand48_set_thread_state(w->conf.par_rng_cur, NULL);

        pthread_mutex_lock(&par->lock);
        if ( ret<0 )
        {
            par->status[ichunk] = -1;
            par->failed = 1;
        }
        else
        {
            par->status[ichunk] = 1;
            par->ndraw[ichunk] = w->conf.par_ndraw;
            par->final[ichunk] = exact || !w->conf.par_ndraw ? 1 : 0;
        }
        pthread_cond_broadcast(&par->cond);
        pthread_mutex_unlock(&par->lock);
    }
    return NULL;
}

// Called with the lock held: once all chunks before a chunk are done, the
// generator state at its start is known. Queue the chunk for processing again
// if it drew random numbers from a wrong state.
static void mplp_par_update_rng(mplp_par_t *par)
{
    while ( par->nrng <= par->nchunk )
    {
        int i = par->nrng - 1;      // the state of chunk i is known
        if ( par->status[i] != 1 ) break;
        if ( !par->final[i] && par->final[i] != -1 )
        {
            par->final[i] = -1;     // queued
            par->rerun[par->nrerun++] = i;
        }
        if ( par->nrng < par->nchunk )
        {
            memcpy(par->rng[par->nrng], par->rng[i], sizeof(par->rng[i]));
            hts_rand48_skip(par->rng[par->nrng], par->ndraw[i]);
        }
        par->nrng++;
    }
}

static void mplp_par_write(mplp_conf_t *conf, mplp_chunk_res_t *res)
{
    htsFile *fp = conf->bcf_fp;
    if ( conf->par_is_vcf==0 )
    {
        if ( bgzf_write(fp->fp.bgzf, res->buf.s, res->buf.l) != res->buf.l ) goto fail;
    }
    else if ( conf->par_is_vcf==1 )
    {
        if ( hwrite(fp->fp.hfile, res->buf.s, res->buf.l) != res->buf.l ) goto fail;
    }
    else
    {
        // vcf_write() flushes BGZF blocks at record boundaries, do the same
        // to produce identical output
        int i;
        char *ptr = res->buf.s;
        for (i=0; i<res->nlen; i++)
        {
            if ( bgzf_flush_try(fp->fp.bgzf, res->len[i]) < 0 ) goto fail;
            if ( bgzf_write(fp->fp.bgzf, ptr, res->len[i]) != res->len[i] ) goto fail;
            ptr += res->len[i];
        }
    }
    return;

fail:
    error("[%s] Error: failed to write to %s\n", __func__,conf->output_fname?conf->output_fname:"standard output");
}

static int mpileup_parallel(mplp_conf_t *conf, bam_hdr_t *hdr, int nregs)
{
    mplp_par_t par;
    memset(&par, 0, sizeof(par));
    par.conf = conf;
    par.hdr = hdr;
    if ( !conf->par_chunk_size ) conf->par_chunk_size = MPLP_PAR_CHUNK_SIZE;

    // Initialize lazily computed data before the threads are started
    uint8_t dummy_seq = 0, dummy_qual = 30, dummy_q = 0;
    int dummy_state = 0;
    probaln_par_t dummy_par = { 0.001, 0.1, 10 };
    probaln_glocal(&dummy_seq, 1, &dummy_seq, 1, &dummy_qual, &dummy_par, &dummy_state, &dummy_q);
    if ( conf->bed )
    {
        int i, nseq;
        char **seqs = regidx_seq_names(conf->bed, &nseq);
        for (i=0; i<nseq; i++) regidx_overlap(conf->bed, seqs[i], 0, 0, NULL);
    }

    mplp_par_plan(&par, nregs);
    if ( !par.nchunk ) { free(par.chunk); return 0; }

    int nworkers = conf->n_threads, i;
    if ( nworkers > par.nchunk ) nworkers = par.nchunk;
    par.max_inflight = 4*nworkers;     // bounds the memory held by finished chunks
    par.res = (mplp_chunk_res_t*) calloc(par.nchunk, sizeof(*par.res));
    par.status = (int*) calloc(par.nchunk, sizeof(*par.status));
    par.final  = (int*) calloc(par.nchunk, sizeof(*par.final));
    par.ndraw  = (uint64_t*) calloc(par.nchunk, sizeof(*par.ndraw));
    par.rerun  = (int*) calloc(par.nchunk, sizeof(*par.rerun));
    par.rng    = (unsigned short (*)[3]) calloc(par.nchunk, sizeof(*par.rng));
    if ( hts_drand48_get_state(par.rng[0]) < 0 ) error("[%s] Error: cannot get the random generator state\n", __func__);
    par.nrng = 1;
    pthread_mutex_init(&par.lock, NULL);
    pthread_cond_init(&par.cond, NULL);

    mplp_worker_t *workers = (mplp_worker_t*) calloc(nworkers, sizeof(*workers));
    for (i=0; i<nworkers; i++) mplp_par_worker_init(&workers[i], &par);
    for (i=0; i<nworkers; i++)
        if ( pthread_create(&workers[i].thread, NULL, mplp_par_worker, &workers[i]) )
            error("[%s] Error: failed to create a thread\n", __func__);

    int ret = 0;
    for (i=0; i<par.nchunk; i++)
    {
        pthread_mutex_lock(&par.lock);
        while ( 1 )
        {
            mplp_par_update_rng(&par);
            if ( par.nrerun ) pthread_cond_broadcast(&par.cond);
            if ( par.status[i] < 0 || par.failed ) break;
            if ( par.status[i] && par.final[i]==1 ) break;
            pthread_cond_wait(&par.cond, &par.lock);
        }
        int status = par.failed ? -1 : par.status[i];
        pthread_mutex_unlock(&par.lock);
        if ( status < 0 ) { ret = -1; break; }

        mplp_par_write(conf, &par.res[i]);
        if ( conf->progress )
        {
            mplp_chunk_t *chunk = &par.chunk[i];
            hts_pos_t pos = chunk->end < hdr->target_len[chunk->tid] ? chunk->end : hdr->target_len[chunk->tid];
            progress_update(conf->progress, hdr->target_name[chunk->tid], pos+1);
        }

        pthread_mutex_lock(&par.lock);
        hts_expand(kstring_t, par.npool+1, par.mpool, par.pool);
        par.pool[par.npool++] = par.res[i].buf;
        memset(&par.res[i].buf, 0, sizeof(par.res[i].buf));
        free(par.res[i].len);
        par.res[i].len = NULL;
        par.nwritten = i+1;
        pthread_cond_broadcast(&par.cond);
        pthread_mutex_unlock(&par.lock);
    }
    pthread_mutex_lock(&par.lock);
    if ( ret < 0 ) par.failed = 1;
    par.finished = 1;
    pthread_cond_broadcast(&par.cond);
    pthread_mutex_unlock(&par.lock);

    for (i=0; i<nworkers; i++) pthread_join(workers[i].thread, NULL);
    for (i=0; i<nworkers; i++) mplp_par_worker_destroy(&workers[i]);
    free(workers);
    for (i=0; i<par.nchunk; i++) { free(par.res[i].buf.s); free(par.res[i].len); }
    for (i=0; i<par.npool; i++) free(par.pool[i].s);
    free(par.pool);
    free(par.res);
    free(par.status);
    free(par.final);
    free(par.ndraw);
    free(par.rerun);
    free(par.rng);
    free(par.chunk);
    pthread_mutex_destroy(&par.lock);
    pthread_cond_destroy(&par.cond);
    return ret;
}
