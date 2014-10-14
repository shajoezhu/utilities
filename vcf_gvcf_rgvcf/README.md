Gernerate .vcf and .seg data from sim-1 population structure.
```bash
../pop_struct.py sim-1 2 1
mv sim-1Samples2msdata1.seg sim-1.seg
```
Use vcf2seg.py to convert .vcf format to .seg format, and compare the result. There should be no difference between the two files.
```bash
./vcf2seg.py -i sim-1Samples2msdata1.vcf -seqlen 30000000
diff sim-1Samples2msdata1.seg sim-1.seg
```

```bash
./vcf2seg.py -i sim-1Samples2msdata1.gvcf -seqlen 30000000
mv sim-1Samples2msdata1.seg sim-1gvcf.seg
```
