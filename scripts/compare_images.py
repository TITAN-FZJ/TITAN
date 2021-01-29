# #!/usr/bin/env python
##################################################################
# Script to compare images
# Adapted from: https://gist.github.com/astanin/626356
##################################################################
import sys   # System library (to read arguments from command line)
from matplotlib.pyplot import imread
import numpy as np

def main():
  file1, file2 = sys.argv[1:1+2]
  # read images as 2D arrays (convert to grayscale for simplicity)
  img1 = to_grayscale(imread(file1).astype(float))
  img2 = to_grayscale(imread(file2).astype(float))
  # compare
  n_m, n_0 = compare_images(img1, img2)
  print("Manhattan norm:", n_m, "/ per pixel:", n_m/img1.size)
  print("Zero norm:", n_0, "/ per pixel:", n_0*1.0/img1.size)
  if(n_m<10000.0) and (n_0<1000.0):
    print(f"Images {file1} and {file2} are equal.")
  else:
    print(f"Images {file1} and {file2} are different!")
    exit(1)


def compare_images(img1, img2):
  # normalize to compensate for exposure difference
  img1 = normalize(img1)
  img2 = normalize(img2)
  # calculate the difference and its norms
  diff = img1 - img2  # elementwise for scipy arrays
  m_norm = np.sum(abs(diff))  # Manhattan norm
  z_norm = np.linalg.norm(diff.ravel(), 0)  # Zero norm
  return (m_norm, z_norm)

def to_grayscale(arr):
  "If arr is a color image (3D array), convert it to grayscale (2D array)."
  if len(arr.shape) == 3:
    return np.average(arr, -1)  # average over the last axis (color channels)
  else:
    return arr

def normalize(arr):
  rng = arr.max()-arr.min()
  amin = arr.min()
  return (arr-amin)*255/rng

if __name__ == "__main__":
  main()


