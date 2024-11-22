#!/opt/homebrew/bin/python3.11

"""
Converts a jpeg/png image to grayscale in pgm format.
Author: Emmanuel A. Larralde Ortiz
"""

import os
import argparse

import cv2


def to_pgm(fname: str) -> None:
    if not os.path.exists(fname):
        return
    pgm_path = f"{os.path.splitext(fname)[0]}.pgm"
    img = cv2.imread(fname)
    gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    cv2.imwrite(pgm_path, gray_img)


def main(args: argparse.ArgumentParser):
    for fname in args.img_files:
        to_pgm(fname)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='Img to PGM',
        description='Converts any common format to pgm'
    )
    parser.add_argument('img_files', nargs='+')

    args = parser.parse_args()
    main(args)

#TESING:
# # any path, even if is not an image, from every workdir.