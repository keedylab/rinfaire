import MainFunctions

def main ():
    multiFlag = True
    args = MainFunctions.setupArguments(multiFlag)
    fileList = MainFunctions.readFile(multiFlag, args)
    networkList = MainFunctions.generateIndividualNetworks(fileList, args)
    multi = MainFunctions.generateMultiNetwork(networkList, args)

if __name__ == "__main__":
    main()