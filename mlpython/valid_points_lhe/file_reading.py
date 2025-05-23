#==============================================================================
# First method looking at each line of the file
#==============================================================================

#
#file_to_open = open("130Demo_out.lha","r") #Open the file
#
#for line in file_to_open:       # Loop through each line
#    line_list = line.split()      # each line is a single string, so we split it into a list of each column or word in that line
#    
#    if line_list[0] == "DECAY" and line_list[1] == "35":
#        print line_list
#        # Or do whatever you want with this line
#
#file_to_open.close()  # Good practice to always close the file when you're done with it


#==============================================================================
# If the file isn't too big, you can read the whole thing into memory as a list
#  and then do your searching through that list. This lets you more easily 
#  access all the data in whatever order you want
#==============================================================================


files_name1 = ["80.lha",
              "90.lha"]
    
files_name2 = ["100.lha",
               "110.lha",
               "120.lha"]


files_name3 = ["130.lha",
              "140.lha",
              "150.lha",
              "160.lha",
              "170.lha",
              "180.lha",
              "190.lha",
              "200.lha",
              "210.lha",
              "220.lha",
              "230.lha",
               "240.lha"]

files_name4 =["250.lha",
              "260.lha",
              "270.lha",
                "280.lha",
               "290.lha",
               "300.lha"]

files_name = files_name4

for name in files_name:
    
    int = files_name.index(name)
    file_to_open = open(name,"r")
    list_of_data_from_file = [line.split() for line in file_to_open]
    file_to_open.close()

    file_to_write = open("Hbb","a")
    
    for i, row in enumerate(list_of_data_from_file):
        if row[0] == "DECAY" and row[1] == "35":
            print(list_of_data_from_file[i+4][0]) #As an example
            file_to_write.write(str(int*10+250))
            file_to_write.write(" ")
            file_to_write.write(list_of_data_from_file[i+4][0])
            file_to_write.write(" ")
            file_to_write.write(list_of_data_from_file[i+4][-2])
            file_to_write.write(" ")
            file_to_write.write(list_of_data_from_file[i+4][-1])
            file_to_write.write("\n")

    file_to_write.close()
