###########################################################################
            if ( not len(exons) ):
                print "   ", ct, tot
###########################################################################


        ct = 0
        tot = 0
        for gene_id in gene_ids:
            tot += 1
            # find all exons associated with the gene id
            exons = find_exons (cursor, gene_id, species)
            if (not exons):
                print  gene_id, " no exons found" 
                exit(1)

            length = 0
            for exon in exons:
                if (not exon.is_constitutive):
                    continue
                length += exon.end_in_gene - exon.start_in_gene + 1

                
            if (length%3):
                ct +=1
                print "\t", gene2stable(cursor, gene_id=gene_id),
                print " number of exons: ", len(exons), 
                print " length ", length, length%3,
                print "   ", ct, tot
                
            # store to gene2exon table
        print " >>> ", ct, "  ", len(exons)



  transcript_id | gene_id | analysis_id | seq_region_id | seq_region_start | seq_region_end 
|           215 |     196 |          57 |            21 |         14769571 |       14800205 |                -1 |
