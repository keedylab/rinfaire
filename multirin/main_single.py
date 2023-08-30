import MainFunctions

def main ():
    multiFlag = False
    args = MainFunctions.setupArguments(multiFlag)
    fileList = MainFunctions.readFile(multiFlag, args)
    networkList = MainFunctions.generateIndividualNetworks(fileList, args)

if __name__ == "__main__":
    main()