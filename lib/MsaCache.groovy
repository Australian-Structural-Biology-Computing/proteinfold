class MsaCache {
    static String directory(Map params) {
        if (params.use_msa_cache != null && !params.use_msa_cache.toString().toBoolean()) return null
        def explicit = params.msa_cache_dir ?: params.boltz_msa_cache_dir ?: System.getenv('PROTEINFOLD_MSA_CACHE_DIR')
        if (explicit) {
            if (!new File(explicit.toString()).isAbsolute()) throw new IllegalArgumentException('msa_cache_dir must be an absolute local path')
            return explicit.toString()
        }
        if (!params.outdir) return null
        if (params.outdir.toString() ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/) return null
        return new File(params.outdir.toString(), 'msa_cache').absolutePath
    }
}
