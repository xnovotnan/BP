# BioVizAnalyzer

## Code Structure and Modularity

BioVizAnalyzer is a dashboard application for viewing and analyzing biomedical data, implemented using Shiny and the R programming language. The application consists of five modules, with each module responsible for a separate data workflow, as shown in the image below:

![Workflows in the interactive dashboard](/docs/Workflow.png)

- [VCFModule](VCFModule/) – analysis of VCF files 
- [VCFComparisonModule](/VCFComparisonModule/) – comparison of multiple VCF files
- [QualimapModule](QualimapModule/) – analysis based on Qualimap outputs
- [QualimapComparisonModule](/QualimapComparisonModule/) – comparison of multiple Qualimap outputs
- [FastqcModule](/FastqcModule/) – analysis of FastQC files

## Setup and Installation

The application can be run in two ways: in [RStudio](#rstudio) or in a [Docker container](#docker-container).

<div id="rstudio">

### Running in RStudio

To run the application, you need to install [all the required libraries](#libraries-version) using RStudio and execute the application with the [app.R file](/app.R) as the entry point.

<div id="libraries-version">

### Libraries Version

| Library      | Version  |
| ------------ | -------- |
| bslib        | 0.9.0    |
| circlize     | 0.4.16   |
| ggridges     | 0.5.6    |
| magrittr     | 2.0.3    |
| patchwork    | 1.3.0    |
| plotly       | 4.10.4   |
| shiny        | 1.10.0   |
| shinyFiles   | 0.9.3    |
| shinyWidgets | 0.9.0    |
| tidyverse    | 2.0.0    |
| tinytex      | 0.56.0   |
| scales       | 1.3.0    |
| fmsb         | 0.7.6    |
| ggrepel      | 0.9.6    |
| png          | 0.1-8    |
| hexbin       | 1.28.5   |

</div>
</div>

<div id="docker-container">
### Running in a Docker Container

Running the application in Docker is easy to set up, as it does not require installing any libraries or having an R compiler installed. To run the app in Docker, execute the `docker-compose.yml` script, which builds the Docker image with the required dependencies and deploys the app locally:

```
docker-compose up --build
```

After a successful build, the application will be available at [localhost:3838](http://localhost:3838).

</div>


## Usage Instructions

The user must open the application in a browser. After opening the application, the user can analyze VCF, Qualimap, or FastQC files. The top panel allows switching between different workflows. Each module contains an upload section where users can select files or folders, an analysis section represented by charts or text fields, and a download section where users can download all the charts from the module. Here is an example of a chart with loaded data:

![Loaded VCF File](/docs/LoadedFile.png)

## Data Privacy and Security

The data in the application is processed locally using files uploaded from the local filesystem. There is no network communication between the application and the internet. The application relies solely on the installed libraries to display charts and can be run entirely offline.
