# Differential Expression Explorer Shiny App for _Amblyomma americanum_

## Quick Start

1. Navigate to [aa-diffex.arcadiascience.com](https://aa-diffex.arcadiascience.com).
2. Log in with your GitHub credentials and authorize the app.
3. Select a model from the upper-most drop-down menu & explore the graphs and analyses available in each tab. See the [Feature of the App](#features-of-the-app) section below for more information.

### Known bugs

* If you refresh the page, it may grey out. Re-navigate to the URL [aa-diffex.arcadiascience.com](https://aa-diffex.arcadiascience.com) and continue exploring.
* GitHub authentication is based on your organization membership. **You must be a member of the Arcadia-Science GitHub organization to get access to the app.** Ping the software team if you are not a member and they can grant you membership.

## Running and Developing Locally

Below we provide instructions for how to run the Shiny app locally.
This is useful if you want to add new functionality to the app.

First, clone the repository to your local computer.
Because this is a private repo, GitHub will authenticate your credentials before it allows you to clone.
If you [followed this lesson](https://training.arcadiascience.com/workshops/20220920-intro-to-git-and-github/lesson/) to set up git on your computer, the authentication should happen automatically.

```
git clone git@github.com:Arcadia-Science/2023-amblyomma-americanum-diffex.git
```

Once you have a copy of the repository, change directories (`cd`) into the `shiny` folder in the repo:

```
cd 2023-amblyomma-americanum-diffex/shiny
```

The Shiny app relies on a few R packages.
To simplify installation of these dependencies, we created an [`environment.yml`](./environment.yml) file that can be used with conda (or mamba) to create an environment for the Shiny app to run in.
If you don't already have conda installed, you can follow the instructions in [this lesson](https://training.arcadiascience.com/arcadia-users-group/20221017-conda/lesson/) to install it.
Once the environment is created, activate it.

```
mamba env create -n diffexshiny -f environment.yml
conda activate diffexshiny
```

The Shiny app relies, on a couple environment variables.
To get started, make a copy of the [.Renviron.dev](/shiny/.Renvrion.dev) and fill in the values.
If you don't have access to the correct values, ping the software team.

```
cd shiny/ && cp .Renviron.dev .Renviron
```

To run the Shiny app, execute the following command:

```
R -e "shiny::runApp('app.R', port = 8100)"
```

This will run the Shiny app from your Terminal.
After the Shiny app is running (~5 seconds), you should see a message that says,

```
Listening on http://0.0.0.0:8100
```

Copy and paste the URL into your browser of choice and the Shiny app should launch.
Even when running locally, the app expects you to authenticate.
When the app launches, you should see a GitHub login screen.
Log in and authorize the app.
You should be redirected to the Shiny app after logging in.

If you have any trouble logging in, try and run the app from a new Private Browser/Incognito window.
Close all other private browsing tabs and windows before opening a new one and copy and pasting the Shiny app URL into the tab.

These instructions assume that you will run the app from the command line.
**Do not try and run the app from RStudio.**
The GitHub authentication doesn't play well with RStudio -- RStudio tries to launch multiple (~100) new tabs for the app when it doesn't immediately launch because of the authentication page.

## Running from a Docker container

First, clone the repository to your local computer.
Because this is a private repo, GitHub will authenticate your credentials before it allows you to clone.
If you [followed this lesson](https://training.arcadiascience.com/workshops/20220920-intro-to-git-and-github/lesson/) to set up git on your computer, the authentication should happen automatically.

```
git clone git@github.com:Arcadia-Science/2023-amblyomma-americanum-diffex.git
```

Once you have a copy of the repository, change directories (`cd`) into the `shiny` folder in the repo:

```
cd 2023-amblyomma-americanum-diffex/shiny
```

Launch Docker Desktop and then navigate to the Terminal and build the Docker container:

```
docker build --progress=plain --platform linux/x86_64 -t diffexshiny .
```

(note that this GitHub repo contains private data, so we have not uploaded the Docker container to Docker Hub.)

Once the Docker container is done building, launch it in a headless state to the 8100 port.

```
docker run -d --platform linux/x86_64 -p 8100:8100 diffexshiny
```

Lastly, navigate to Docker Desktop and click the "Containers" tab.
Launch the app by clicking `6626:6626` in the "Port(s)" column.

## Overview of the App

This Shiny application serves as an interactive tool for exploring differential gene expression results derived from two differential expression models built from publicly available short-read RNA-seq data for the lone star tick _Amblyomma americanum_.

## Data Preparation

The Shiny app is built from data files that are packaged in the [`input_data`](./input_data) folder in this repository.
The folder contains metadata associated with the experiment, DESeq2 differential expression models (see the [Snakefile](../Snakefile) in this repository), and [ortholog and functional annotations for the _A. americanum_ genome](https://github.com/Arcadia-Science/protein-data-curation).

## Features of the App

### Model Selection

Users can select one of two differential expression models to explore: `sex_tissue` or `sex_tissue_blood_meal_hour`.
The descriptions and the rationale behind these models are provided within the application.

### PCA Plot

This tab provides a Principal Component Analysis (PCA) plot based on variance stabilized transformed data, where users can color the points based on selected metadata variables like sex, tissue, blood meal hour range, or study title.

### Differential Expression (DE) Analysis

In this tab, users can filter the differential expression results based on log2 Fold Change, p-value, and base mean values.
They can select conditions for comparison, and optionally upload a CSV file to include additional information for coloring the points in the output plot.
The results can be visualized in a Volcano and MA plot.
Additionally, users can download the filtered results as a CSV file.

### Expression of Specific Genes

Here, users can enter a specific gene name to generate a boxplot showing its expression across different conditions.
Each boxplot is overlaid with the scatter points representing the expression of individual samples.

### Expression by Condition

This tab reports on expression per condition, allowing users to filter genes based on their expression thresholds and percentiles.
Users can view expression distribution and download tables of genes that are never, sometimes, or always expressed in the selected condition.

## Troubleshooting and issues

Sometimes the Shiny app may freeze.
In testing, this happened most frequently on the DE Analysis tab when trying to switch which data was plotted (e.g. switching models or adding user data).
If this happens:

* **When running in the browser:** exit the tab and re-launch the URL in a new tab.
* **When running locally:** the best thing to do is to exit the browser tab, kill the Shiny app by pressing <kbd>ctrl</kbd> + <kbd>c</kbd>, and to relaunch it with the R command `R -e "shiny::runApp('app.R', port = 8100)"`.

Then, make sure the first data you plot in the DE Analysis tab is the data that froze the app.

## Next Steps

We don't have any next steps planned for the app.
If you have ideas about missing functionality or data, talk to the software team or file and issue on the repo.
