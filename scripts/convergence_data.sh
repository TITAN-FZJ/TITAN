# This file is meant to be used to extract the convergence history
# of the self consistency even while TITAN is running
# the suggested use is as follows
#
# $./convergence_history.sh path_to_file >> output_file
# 
# in this case it extracts the data and writes it in a file named "output_file"
#
# If you only want to print the values in the terminal then omit the last part, as in
# 
# $./convergence_history.sh path_to_file
#

# Gets the path to the file
file="$1"

# fvec is an array with several components, I am only interested on the last one
# the command below does the following
# - looks for all the ocurrences of the word fvec
# - takes only the last line (this line will contain always the last component of fvec)
# - gets the text fvec( x) and stores it in the variable "tracker"
# fvec has a different number of components on different simulation, but this method
# makes sure you always get the text correpondig to the last component
tracker=`grep -i "fvec( " output/selfconsistency | tail -n 1 | cut -c 1-26`

# looks only for the lines where we print the last component of fvec
data=`grep -i "$tracker" $file`

# removes the "fvec( x)" text from the data
cleanData=`printf "$data\n" | sed "s/$tracker//g"`

# prints only the number associated with the last component of fvec
echo $cleanData
