import re

#initialize required variables
newline = 0 #will check if a newline needs to be added to the string which will occur at every 60th nucleotide
line = 0 #will vary between 0 and 1 to indicate false or true and will send this value to the add space function
times_reached = 0 #counts the number of times that a newline has been added in the given for loop
count = 10 # counts the index at which the space or newline should be added
previous_splice = 0 #the previous splice position
add = 1 #add previous nucleotide index to new index

#allows the user to input the name of the fasta file
fasta_name = input("please input the name of the fasta file : ")


#opens the file whos name was provided by the user
fasta_file = open(fasta_name)

#reads the content of fasta_file line by line
fasta_lines = fasta_file.readlines()

#saves the header of the fasta sequence
fasta_header = fasta_lines[0]

#deletes the header from the fasta_file, then removes all special characters and saves the fasta sequence as a single line
fasta_lines.pop(0)
fasta_sequence = ''.join(fasta_lines)
re.sub(r"[\n\t\s1-9]*", "", fasta_sequence)

#allows the user to input the name of enzyme text file
enzyme_file_name = input("please input the name of the enzyme text file : ") + '.txt'

#opens the file whos name was provided by the user
enzyme_file = open(enzyme_file_name)

#reads the content of the enzyme text file line by line 
enzyme_list = enzyme_file.readlines()

#this for loop will separate the enzyme name from the location of its splice site and save each of them to a list
for x in range(len(enzyme_list)):
	#if this is the first loop, then lists must be initialized
	if x == 0:
		#the enzyme list will be split into to columns by the ; character and saved to two different lists
		enzyme_name = [enzyme_list[0].split(";")[0]]
		enzyme_site = [enzyme_list[0].split(";")[1]]
		enzyme_site[0] = ''.join(enzyme_site[0].split())
	#if this is not the first loop, append the previous lists
	else:
		enzyme_name.append(enzyme_list[x].split(";")[0])
		enzyme_site.append(enzyme_list[x].split(";")[1])
		enzyme_site[x] = ''.join(enzyme_site[x].split()) 

#create an add space function which will be called to either add a space to the nucleotide character or to add a newline if the 60th character has been reached
def add_space (string, integer, line,multiply, add):
	#line checks if a space or a newline is to be added (line = 0 means space and line = 1 means a newline)
	if line == 0:
		#return a string which has a space added at the index given by the integer variable
		return string[:integer] + ' ' + string[integer:]
	#if line = 1
	else:
		#adds a newline at the integer index, then prints the nucleotide position and adds a tab character before continuing the nucleotide string
		return string[:integer] + '\n' +str(round(60*multiply+add)) + '\t' + string[integer:]

#This for loop will be used to create the lists which contain the nucleotide strings with their nucleotide position as well as the spaces and newlines added
for x in range(len(enzyme_name)):
	#checks if the enzyme splice site is present in the fasta sequence before commiting to all the other loops
	if fasta_sequence.find(enzyme_site[x].replace('%', '')) != -1:
		#determines the fasta_splice_position by finding the indeces in the fasta file at which the enzyme splice sequence is found 
		fasta_splice_position = [m.start() for m in re.finditer(enzyme_site[x].replace('%', ''), fasta_sequence)]
		#if this is the first instance of the loop, the splice_counter list needs to be initialized with the number of splice positions found
		if x == 0:
			splice_counter = [len(fasta_splice_position)]
		#otherwise, append the splice_counter list with the number of splice sites for the given enzyme
		else:
			splice_counter.append(len(fasta_splice_position))
		#locates the index of the % (representing the splice site) in the given enzyme
		splice_enzyme = enzyme_site[x].find("%")
		#This loop will split the nucleotide string based on the splice position
		for i in range(0, len(fasta_splice_position)):
			#The index at which the space or newline will be added 
			split_at = fasta_splice_position[i] + splice_enzyme - previous_splice
			#if this is the first loop, lists must be initialized
			if i == 0:
				#string containing the right side of the fasta sequence at the splice location
				spliced = fasta_sequence[:split_at]
			else:
				#equal to the right side of the index location of the new fasta string
				spliced = new_fasta[:split_at]
				#equal to the left side of the index location of the previous fasta string
				new_fasta = new_fasta[split_at:]
				#this will add the value of the previous splice index 
				add = fasta_splice_position[i-1]+splice_enzyme+1
			#This loop will add spaces and newlines where they are required based on the index value
			for j in range(int(len(spliced)/10)):
				#newline will increment until 6 at which point a newline will be added
				newline = newline + 1
				if newline == 6:
					#resets newline
					newline = 0
					line = 1
					#sums the number of times that this nucleotide segment has added a newline 
					times_reached = times_reached + 1
				#calls the add space function
				spliced = add_space(spliced, count, line, times_reached, add)
				#count determines the index location at which the space will be added
				if line == 1:   count = count + len(str(times_reached*60+add))+1
				count = count + 11
				line = 0
			#reset all values for next loop
			newline = 0
			times_reached = 0
			count = 10
			#if this is the first instance of the loop, the spliced fasta list and new fasta string are initialized
			if i == 0:
				spliced_fasta = [spliced]
				new_fasta = fasta_sequence[split_at:]
			else:
				spliced = str(fasta_splice_position[i-1]+splice_enzyme+1) + '\t' + spliced
				spliced_fasta.append(spliced)
			#previous splice contains the previous index at which point the string was spliced
			previous_splice = split_at
		#Performs the loop to add  spaces and newlines one more time on the right-most segment of the  spliced fasta sequence
		for j in range(int(len(new_fasta)/10)):
			newline = newline + 1
			add = fasta_splice_position[len(fasta_splice_position)-1]+splice_enzyme+1  
			if newline == 6:
				newline = 0
				line = 1
				times_reached = times_reached + 1
			new_fasta = add_space(new_fasta, count, line, times_reached, add)
			if line == 1:   count = count + len(str(times_reached*60 + add))+1
			count = count + 11
			line = 0
		#resets all values for next loop
		newline = 0
		count = 10
		times_reached = 0
		new_fasta = str(fasta_splice_position[len(fasta_splice_position)-1]+splice_enzyme+1) + '\t' + new_fasta
		spliced_fasta.append(new_fasta)
		#all_splices contains lists of all segments from all enzymes
		if x == 0:
			all_splices = [spliced_fasta]
		else:
			all_splices.append(spliced_fasta)
	#if the fasta sequence did not contain the sequence that is spliced by the enzyme
	else:
		#initializes the lists but sets their first value to none or 0
		if x == 0:
			splice_counter = ["none"]
			all_splices = ["0"]
		else:
			splice_counter.append("none")
			all_splices.append("0")
	#resets values for next loop
	previous_splice = 0
	add = 1

#prints all results
print ("Restriction enzyme analysis of sequence from file " + fasta_name)
print ("Cutting with enzymes found in the file " + enzyme_file_name)
print ("----------------------------------------------------------------------------------------------")
print ("sequence name: " + fasta_header.replace('>', '').replace('\n', ''))
print ("Sequence length is " + str(len(fasta_sequence)) + " bases long")
print ("----------------------------------------------------------------------------------------------")
for x in range(len(enzyme_list)):
	if str(splice_counter[x]) == "none":
		print ("There are no cutting sites for " + str(enzyme_name[x]))
		print ("\n---------------------------------------------------------------------------")
	else:
		print ("There are " + str(splice_counter[x]) + " cutting sites for " + str(enzyme_name[x]) + ", cutting at " + str(enzyme_site[x]))
		print ("There are " + str(splice_counter[x]+1) + " fragments\n")
		for i in range(splice_counter[x]+1):
			print ("Length: " + str(len(re.sub(r"[\n\t\s1-9]*", "", all_splices[x][i]))))
			if i == 0:
				print ('1\t' + str(all_splices[x][i]))
			else:
				print (str(all_splices[x][i]))
		print ("\n\n----------------------------------------------------------------------------------------------")

