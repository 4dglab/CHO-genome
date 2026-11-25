# PURPOSE: nucmer coords -> chromosome assignment
# INPUT:   coords file
# OUTPUT:  chr_assignment.tsv
# TOOL:    awk

BEGIN {
    FS = "\t"
}
NR <= 4 { next }   # skip header
{
    ref   = $12     # reference (chr)
    qry   = $13     # query (scaffold)
    alen  = $8      # alignment length

    # accumulate alignment bases per (qry, ref) pair
    key = qry SUBSEP ref
    if (!(key in bases)) {
        bases[key] = 0
    }
    bases[key] += alen

    # record all scaffold names
    scaffolds[qry] = 1
}
END {
    # for each scaffold, find the chromosome with max total aligned bases
    for (s in scaffolds) {
        max_bases = 0
        best_chr  = "unassigned"
        total_aligned = 0

        for (r in chrs) {
            key = s SUBSEP r
            b = bases[key]
            if (b > 0) {
                total_aligned += b
                if (b > max_bases) {
                    max_bases = b
                    best_chr  = r
                }
            }
        }
        printf "%s\t%s\t%d\n", s, best_chr, total_aligned
    }
}
