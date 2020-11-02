#!/usr/bin/env python
import sys   # System library (to read arguments from command line)

# Defining the tolerance
tol = 1.e-7

################################################################################
# Checks if a string is a number
################################################################################
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

################################################################################
# Get the data from the file and save it into a list
################################################################################
def read_data(filename):
  # open a file using with statement
  data = []
  with open(filename,'r') as file:
    for line in file:
      # check if the current line starts with "#" or "!"
      if not (line.startswith("#") or line.startswith(" #") or line.startswith("!")):
        data_temp = [float(item) for sublist in [x.split() for x in line.split("=")] for item in sublist if is_number(item)]
        data.append(data_temp)
  return len(data),data


if __name__ == "__main__":
  # Getting filename from arguments of command line
  file1 = sys.argv[1]
  # Reading file 1
  nlines1, data1 = read_data(file1)
  # Getting filename from arguments of command line
  file2 = sys.argv[2]
  # Reading file 2
  nlines2, data2 = read_data(file2)

  # Getting tolerance from command line, when given
  if(len(sys.argv)>3):
    tol = float(sys.argv[3])

  # Check if number of lines is the same
  if ( nlines1 != nlines2 ):
    print("Number of lines in file 1: {}".format(nlines1))
    print("Number of lines in file 2: {}".format(nlines2))
    print("Files have different number of lines")
    exit(1)

  # Loop over the lines of the files (that are equal due to the check above)
  err  = 0 # Start error counter
  line = 1
  for line1,line2 in zip(data1,data2):
    # Number of elements in line i of file 1
    ncols1=len(line1)
    # Number of elements in line i of file 2
    ncols2=len(line2)

    if ( ncols1 != ncols2 ):
      print("Number of columns in line {} of file 1: {}".format(line,ncols1))
      print("Number of columns in line {} of file 2: {}".format(line,ncols2))
      print("Files have different number of columns")
      exit(1)

    # Loop over the columns of line i (that are equal due to the check above)
    col = 1
    for value1,value2 in zip(line1,line2):
      if ( abs(value1 - value2) > tol ):
        print("Column {} of line {} has different values: {} {}, difference: {}".format(col,line,value1,value2,abs(value1-value2)))
        err = err + 1
      col = col + 1

    line = line + 1

if(err == 0):
  print("All values are equal inside the tolerance {}!".format(tol))
  exit(0)
else:
  print("Different values found greater than tolerance {}!".format(tol))
  exit(1)


