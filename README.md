
### BEF_rock_pools

All code and data to reproduce the analysis reported in the following publication which is currently in preparation can be found in this repository:

+ Snoeks JM, Pathak-Broome A, Hagan JG & Vanschoenwinkel B (in prep.). Global variation in biodiversity–ecosystem functioning (BEF) relationships in aquatic metacommunities on terrestrial islands.

#### Download the repository

To reproduce the analyses reported, first download this repository to your computer. This can be done in two ways:

**1. with git**

For this, in the Terminal (Mac) navigate to the folder where you want the repository to live (commands may differ on Windows)

+ cd path/to/local/folder

Then, type the following command:

git clone https://github.com/JorenMSnoeks/BEF_rock_pools

This should download the repository.

**2. without git**

If you don't have git installed, you can download the repository as zip file and save it locally.

--> Clone or download (green button top right) --> Download Zip

then save and extract the zip where you want the directory to be saved on your computer. To run the code correctly, it is important to create a R-Project. In R-Studio go to File > New Project... > Existing Directory > Choose the extracted directory.

#### Download the relevant packages

For package version control, we used *renv* (Ushey and Wickham 2023). Therefore, in the repository, there is a renv lockfile (renv.lock) that contains information about the versions of all packages used in the analysis. 

To download the correct versions of all packages that we used in the analysis, you need to install the renv package:

+ install.packages("renv")

Once renv is installed, in the R-console, run the following command:

+ renv::restore()

The restore() command will download the correct versions of all packages used in the analysis into a local folder in the repository called 'renv'.

#### Run the scripts

In the 'Code' folder, there are 15 Quarto notebooks. Running these in order will reproduce the analyses presented in the paper and provide additional context on the data and methods.






