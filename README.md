# DBCscreen
A pipeline to rapidly detect DNA barcode contamination for biodiversity research

## To ensure DBCscreen runs smoothly, you should first complete the following three steps:

1. Download [FCS-GX](https://github.com/ncbi/fcs) 
2. Download the [DBCscreen database](https://doi.org/10.57760/sciencedb.20428) 
3. Update the path in the [config file](https://github.com/xiebio/DBCscreen/blob/main/DBCscreen.config) to point to the downloaded database and gx.

## Script for DNA barcode contamination screening

```PYTHON
python DBCscreen.py -i input_file -o ./results
```

Any suggestions or problem, please contact Jiazheng Xie（xiejz@cqupt.edu.cn) .

*Citation*：
to be added
