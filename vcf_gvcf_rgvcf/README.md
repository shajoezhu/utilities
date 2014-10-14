../pop_struct.py sim-1 2 1

mv sim-1Samples2msdata1.seg sim-1.seg

vcf2seg.py -i sim-1Samples2msdata1.vcf -seqlen 30000000
diff 1Samples2msdata1.seg sim-1.seg

``~bash
vcf2seg.py -i sim-1Samples2msdata1.gvcf -seqlen 30000000
diff 1Samples2msdata1.seg sim-1.seg
``
