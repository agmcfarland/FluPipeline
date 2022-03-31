# runner script to programatically create sample and run Rmarkdown reports

library(optparse)

option_list = list(
  make_option(c("--softwareDir"), type="character", default='/home/common/SARS-CoV-2-Philadelphia', help="path to software directory", metavar="character"),
  make_option(c("--report_type"), type="character", default='sample', help="type of report to generate", metavar="character"),
  make_option(c("--sampleDir"), type="character", default='sample directory', help="sample directory path", metavar="character"),
  make_option(c("--samplename"), type="character", default='sample name', help="sample name", metavar="character"),
  make_option(c("--baseDir"), type="character", default='base directory', help="type of report to generate", metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#troubleshooting inputs sample start
# opt$sampleDir <- '/home/agmcfarland/quick_tests/output/sampleOutputs/Ashley_1_2_S81'
# opt$samplename <- basename(opt$sampleDir)
# opt$report_type <- 'sample'
# opt$softwareDir <- '/home/agmcfarland/flu_project/FluPipeline'
#troubleshooting inputs sample end

#troubleshooting inputs run start
# opt$softwareDir <- '/home/agmcfarland/flu_project/FluPipeline'
# opt$report_type <- 'run'
# opt$baseDir <- '/home/agmcfarland/flu_project/FluPipeline/run_test/output'
#troubleshooting inputs run end

if (opt$report_type=='sample'){
  suppressWarnings(rmarkdown::render(file.path(opt$softwareDir, 'sample_report.Rmd'),
                    output_file = file.path(opt$sampleDir, paste0(opt$samplename, '.pdf')),
                    params = list(
                      date  = format(Sys.time(), "%Y-%m-%d"),
                      title = opt$samplename,
                      sampleDir=opt$sampleDir,
                      samplename=opt$samplename),
                    clean=TRUE))
  } else if (opt$report_type=='run') {
    suppressWarnings(rmarkdown::render(file.path(opt$softwareDir, 'run_report.Rmd'),
                    output_file = file.path(opt$baseDir, 'runSummary.pdf'),
                    knit_root_dir=opt$softwareDir,
                    params = list(
                      date  = format(Sys.time(), "%Y-%m-%d"),
                      title = paste0('FluPipeline Run Summary'),
                      baseDir = opt$baseDir),
                    clean=TRUE))
    } else {
  system('echo Error in report_runner.R. Check that report_type specified is correct/exists')
    }





