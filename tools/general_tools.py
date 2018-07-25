from __future__ import print_function

AtomNameToNumber= {'H':  1, 'He': 2,
'Li': 3, 'Be': 4, 'B':  5, 'C':  6,  'N': 7,  'O': 8, 'F':  9, 'Ne':10,
'Na':11, 'Mg':12, 'Al':13, 'Si':14,  'P':15,  'S':16, 'Cl':17, 'Ar':18,
'K': 19, 'Ca':20, 
'Sc':21, 'Ti':22, 'V': 23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cd':29, 'Zn':30,
'Ge':31, 'Ga':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36, 
'Rb':37, 'Sr':38,
'Y':39,  'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48,
'In':49, 'Sn':50, 'Sb':51, 'Te':52,  'I':53, 'Xe':54,
'Cs':55, 'Ba':56,
'La':57, 'Hf':72, 'Ta':73,  'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80,
'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86
}

AtomNumberToName =  dict((value,key) for key,value in AtomNameToNumber.iteritems())


def readFile(fileName,option="rb+"):
    """Read File and return text in file
       
       fileName = str, Name of the file to be read (or path+fileName)
       return   = str, Information in fileName
    """

    try:
        data = open(fileName,option)
    except:
        raise Exception(" Error reading file: " + fileName)

    text = data.read()
    data.close()  
    return text

def writeOutput(fileName,content):
    """ write content to file [fileName]

        fileName = str, Name of the file to be read 
        content  = str, content written to the file
        return   = None
    """
    try:
      OUT = open(fileName,"w")
    except:
      raise Exception("Error writing to file: "+fileName)
    
    OUT.write(content)
    OUT.close()
    return;
def getMoldenAtoms(i,coords):
        return "% 4s %4d  %4d    % 14.10f   % 14.10f   % 14.10f \n"  % (AtomNumberToName[coords[0]],i, coords[0], coords[1],coords[2],coords[3])

def scaleCoords(Coords,scale=1.889725989):
    """ Assumes Coords to be a list of [[ Atom , x , y ,z ],...] entries 
        Default scale = Angstrom to Bohr conversion """
        
    return [ [Coord[0],scale*Coord[1],scale*Coord[2],scale*Coord[3]] for Coord in Coords ]

def GetElementsByLine(text,ElementPositions,Vartype=[str]):
    """ Get Elements from a text by line, it is assumend that in every line looks the same.  
        it splits the text by line=text.splitlines() and reads the elements of ElementPostions of line.split()
       
        text                = str          
        ElementPositions    = int/list(int),The Elements that should be read in every Line of Text
                                            e.g. [2,3,4] reads the 3,4,5 element of the lines in the text
                                            if ElementPositions = -1 the whole line is read, therefore the Vartype needs to  be a string!
        Vartype             = list(type)   ,List of types, either with just one entry, or with as many as ElementPositions has! 
                                            gives the type of the indivual variable that is read by list
                                            e.g. [str,float,float] would mean taht the 3rd element of every line (splited by spaces)
                                            will be read as a string, the 4th and 5th as a float
                                            if Vartype has only one entry, all elements in ElementPositons are read as the type

        return             = list(list(Vartype)), list of elements. first list has Number of Lines elements and every entry has len(ElementPositions) elements 
    """
        
    Elements=[]
    if text is None:
        return
    if ElementPositions == -1 and Vartype[0] != str:
       print("Attention, Vartype = " + str(Vartype) + " but whole line should be read! Changing type to str")
    if type(ElementPositions) == int:
        for line in text.splitlines():
            if ElementPositions !=-1:
                column=line.split()
                Elements.append(Vartype[0](column[ElementPositions]))
            else:
                Elements.append(line)
        return Elements
    elif type(ElementPositions) == list:
        for line in text.splitlines():
           column=line.split()
           elements = []
           if len(Vartype) == 1:
              Vartype = len(ElementPositions)*Vartype
           for i  in  range(len(ElementPositions)):
               if type(ElementPositions[i]) == int:
                  elements.append(Vartype[i](column[ElementPositions[i]]))
               else:
                  raise Exception("ElementPositions can only be a list of integers!!")
           Elements.append(elements)
        return Elements 
    else: 
            raise Exception("ElementPositions can only be integer or a list of integers!!")

def getValues(text,Keywords,Length,Shift,ElementPositions,Vartype=[str]): 
    """ reads Values at postion of ElementPositions py lines from a file
        it uses for this the function NFind() to get the text and 
        GetElementsByLine() 
        
        text                = str          ,Text fromFile that was read, 
                                            NFind(text,Keyword,Length,Shift) has to give a fragment of Text that has
                                            every line the same structure or at least N elements seperated by spaces, 
                                            where N is the maximal element in ElementPositions 
        Keywords            = str/list(str),One or a list of Keywords that are in the string
        Length              = int/list(int),Int or list of ints, gives the length in Characters of the fragment text
        Shift               = int/list(int),Gives the shift in Characters from the keyword everything can be read

        Keywords,Length and Shift need to be of the same number of elements! In the right order! One can use the getParam program
        to get the right parameters, in a very easy way.

        ElementPositions    = int/list(int),The Elements that should be read in every Line of the fragment of the Text
                                            e.g. [2,3,4] reads the 3,4,5 element of the lines in the fragment text
                                            if ElementPositions = -1 the whole line is read, therefore the Vartype needs to  be a string!
        Vartype             = list(type)   ,List of types, either with just one entry, or with as many as ElementPositions has! 
                                            gives the type of the indivual variable that is read by list
                                            e.g. [str,float,float] would mean taht the 3rd element of every line (splited by spaces)
                                            will be read as a string, the 4th and 5th as a float
                                            if Vartype has only one entry, all elements in ElementPositons are read as the type

        return              = list(Vartype),list of elements, where every line is a own list!
    """

    if type(Keywords) == str:
        grep=NFind(text,Keywords,Length,Shift) 
        return GetElementsByLine(grep,ElementPositions,Vartype)
    elif type(Keywords) == list:
        grep=NFind(text,Keywords,Length,Shift) 
        Elements=[]
        for i in range(len(grep)):
            out=GetElementsByLine(grep[i],ElementPositions[i],Vartype)
            Elements += out
        return Elements

def NFind(text,Keywords,Length,Shift): 
    """ NFind calls Nfind() to get multiple greps.  
        Nfind() uses str.find()
        
       
        text                = str          ,Text should have the keyword inside 
        Keywords            = str/list(str),One or a list of Keywords that are in the string
        Length              = int/list(int),Int or list of ints, gives the length in Characters of the fragment text
        Shift               = int/list(int),Gives the shift in Characters from the keyword everything can be read

        Keywords,Length and Shift need to be of the same number of elements! In the right order! One can use the getParam program
        to get the right parameters, in a very easy way.

        return             = str/list(str), part of the text that was greped, with Shift from the Keyword and Lenght Character's
    """
    #
    #  equal to NFind but also possible with multiple keywords
    #  ATTENTION: works only if the KEYWORDS are in the correct order!!
    #  
    #
    output=[] 
    end=0
    if type(Keywords) == str:
       output,end=Nfind(text,Keywords,length=Length,shift=Shift,begin=end)
    elif type(Keywords) == list:
       for i in range(len(Keywords)):
           out,end=Nfind(text,Keywords[i],length=Length[i],shift=Shift[i],begin=end) 
           output.append(out)
    else:
       print("Keywords can only be string or a list of strings")
       sys.exit()
    return output

def Nfind(text,keyword,length=100,shift=0,begin=0):
    """ Nfind() greps part of a text, the fragment Text has lenght characters and is the text in shift characters from the keyword 
        uses the str.find() function of python, begin gives in characters the start where it will start to search for the keyword.
         
        
       
        text                = str          ,Text should have the keyword inside 
        keyword             = str          ,One Keyword that is the string 
        length              = int          ,Int or list of ints, gives the length in Characters of the fragment text
        shift               = int          ,Gives the shift in Characters from the keyword everything can be read

        Keywords,Length and Shift need to be of the same number of elements! In the right order! One can use the getParam program
        to get the right parameters, in a very easy way.

        return             = str           ,part of the text that was greped, with Shift from the Keyword and Lenght Character's
    """
    # 
    # Searches for keyword in text and returns the 
    # len characters (100=default) starting from the
    # beginning of the keyword, alternative can shift the print 
    #
    start=text.find(keyword,begin)
    if start == -1:
       print(keyword + " not found in text!")
       return None,begin
    if start+shift+length > len(text):
       return text[start+shift:len(text)],-1
    else:
       return text[start+shift:start+shift+length],start+shift+length
