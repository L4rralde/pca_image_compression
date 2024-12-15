# Image compression with dominant eigenvectors and interpolation

## Building

```sh
mkdir cpp/output/
!g++ cpp/include/*/*.c* cpp/src/main.cpp -o cpp/output/main.o
```

## Use

```sh
./cpp/output/main.o <src_img_path> <dst_img_path> [<k>]
```

For example,

```sh
./cpp/output/main.o images/HmrZq.pgm result.pgm 128
```
