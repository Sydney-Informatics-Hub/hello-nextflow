# Supporting materials

## Recommended resources

Here are some useful resources we recommend to help you get started with running nf-core pipelines and developing Nextflow pipelines:

### Developed by us

* [SIH Nextflow template](https://github.com/Sydney-Informatics-Hub/template-nf)
* [SIH Nextflow template guide](https://sydney-informatics-hub.github.io/template-nf-guide/)
* [SIH Customising nf-core workshop](https://sydney-informatics-hub.github.io/customising-nfcore-workshop/)
* [Australian BioCommons Seqera Platform Service](https://www.biocommons.org.au/seqera-service )
* [NCI Gadi nf-core instutitonal config](https://nf-co.re/configs/nci_gadi/)
* [Pawsey Setonix nf-core instutitional config](https://nf-co.re/configs/pawsey_setonix/)

### Developed by others 

* [Nextflow training](https://training.nextflow.io/)
* [Nextflow patterns](https://nextflow-io.github.io/patterns/index.html)
* [Nextflow blog](https://www.nextflow.io/blog.html)
* [Nextflow coding best practice recommendations](https://carpentries-incubator.github.io/Pipeline_Training_with_Nextflow/07-Nextflow_Best_Practice/index.html)
* [Seqera community forums](https://community.seqera.io/)

## Nextflow tips and tricks

Nextflow has some useful features for executing pipelines and querying metadata and history. Here are some resources to help you get started.

### Query specific pipeline executions

The [Nextflow log](https://www.nextflow.io/docs/latest/cli.html#log) command is useful for querying execution metadata and history. You can filter your queries and output specific fields in the printed log. 

```default
nextflow log <run_name> -help
```

### Execute Nextflow in the background

The [`-bg`](https://www.nextflow.io/docs/latest/cli.html?highlight=bg#execution-as-a-background-job) options allows you to run your pipeline in the background and continue using your terminal. It is similar to `nohup`. You can redirect all standard output to a log file. 

```default 
nextflow run <workflow_repo/main.nf> -bg > workshop_tip.log
```

### Capture a Nextflow pipeline's configuration

The [Nextflow config](https://www.nextflow.io/docs/latest/cli.html#config) command prints the resolved pipeline configuration. It is especially useful for printing all resolved parameters and profiles Nextflow will use to run a pipeline. 

```default
nextflow config <workflow_repo> -help
```

### Clean Nextflow cache and work directories

The [Nextflow clean](https://www.nextflow.io/docs/latest/cli.html#clean) command will remove files from previous executions stored in the `.nextflow` cache and `work` directories. The `-dry-run` option allows you to preview which files will be deleted. 

```default
nextflow clean <workflow_repo> -help
```

### Change default Nextflow cache strategy
Workflow execution is [sometimes not resumed as expected](https://training.nextflow.io/basic_training/cache_and_resume/#resume-troubleshootingl). The [default behaviour of Nextflow cache keys](https://www.nextflow.io/docs/latest/process.html#cache) is to index the input files meta-data information. Reducing the cache stringency to `lenient` means the files cache keys are based only on filesize and path, and can help to avoid unexpectedly re-running certain processes when `-resume` is in use. 

To apply lenient cache strategy to all of your runs, you could add to a custom configuration file:

```default
process {
    cache = 'lenient'
}
```

You can specify different cache stategies for different processes by using `withName` or `withLabel`. You can specify a particular cache strategy be applied to certain `profiles` within your institutional config, or to apply to all profiles described within that config by placing the above `process` code block outside the `profiles` scope.    

### Access private GitHub repositories

To interact with private repositories on GitHub, you can provide Nextflow with [access to GitHub](https://www.nextflow.io/docs/latest/sharing.html#github-credentials) by specifying your GitHub user name and a [Personal Access Token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) in the [`scm` configuration file](https://www.nextflow.io/docs/latest/sharing.html#scm-configuration-file) inside your specified `.nextflow/` directory:

```default
providers {

  github {
    user = 'georgiesamaha'
    password = 'my-personal-access-token'
  }

}
```
