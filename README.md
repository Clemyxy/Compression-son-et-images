# CTC

Only works on Linux or Windows Subsystem for Linux.


### How to compile

go to the TP_DCT folder, then use the following command:

```bash
make
```

to run the tests.

To create the pages showcasing the results of the compression algorithms use one of the following commands:

##### Linux

```bash
./page_jpeg
```
```bash
./page_son
```
```bash
./page_ondelette
```

Once ran, you will find a file named xxx.html in the TP_DCT folder, open it in you internet browser to see the results.

##### WSL

if you are using wsl, you'll need to first use the following commands:

```bash
sudo apt-get update
sudo apt-get install netpbm
```

Then you can follow linux instructions to generate the page you want.

