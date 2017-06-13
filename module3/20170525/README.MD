
# How to convert CSV files to NWB format
> By Adolfo Hoyos


From the path:


```
module3/20170525
```

Take [this file](main.r) and [this file](to.nwb.r), and download them to the same path.

Make yourself sure that any **file to be converted has its first row containing column names**, and then run:

```
$ ./main.r <filename>
```

A file will be generated in that path with the same name and **.nwb** extension.



	
	