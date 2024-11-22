#!/opt/homebrew/bin/python3.11

import os, sys, random, glob, subprocess
import pytest

sys.path.append('../python')

import img_to_pgm as pgm

def call_img_to_pgm():
    path = os.getcwd()
    root_path = subprocess.check_output(["git", "root"]).decode().strip()
    filename = random.choice(glob.glob(f"{root_path}/*"))
    if os.path.isfile(filename):
        os.chdir(root_path)
    else:
        os.chdir(filename)
    image_file = random.choice(glob.glob(f"{root_path}/*/*"))
    pgm.to_pgm(image_file, test=True)
    os.chdir(path)

@pytest.mark.parametrize('execution_number', range(100))
def test_img_to_pgm(execution_number):
    call_img_to_pgm()
