from generate import MainFunctions
import logging

def main ():
    multiFlag = True
    args = MainFunctions.setupArguments(multiFlag)
    logging.basicConfig(filename=f'{args.output}multirin.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s', level=logging.INFO)
    fileList = MainFunctions.readFile(multiFlag, args)
    networkList = MainFunctions.generateIndividualNetworks(fileList, args)
    multi = MainFunctions.generateMultiNetwork(networkList, args)

if __name__ == "__main__":
    main()