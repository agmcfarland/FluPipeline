# runner script to programatically create sample and run Rmarkdown reports

library(optparse)

option_list = list(
  make_option(c("--softwareDir"), type="character", default='/home/common/SARS-CoV-2-Philadelphia', help="path to software directory", metavar="character"),
  make_option(c("--report_type"), type="character", default='sample', help="type of report to generate", metavar="character"),
  make_option(c("--sampleDir"), type="character", default='sample directory', help="sample directory path", metavar="character"),
  make_option(c("--samplename"), type="character", default='sample directory', help="sample name", metavar="character"),
  make_option(c("--baseDir"), type="character", default='base directory', help="type of report to generate", metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#troubleshooting inputs sample start
#opt$softwareDir <- '/home/agmcfarland/flu_project/flu_pipeline'
#opt$report_type <- 'sample'
#opt$samplename <- 'IBV_Yamagata_Ref_snpindel'
#opt$sampleDir <- '/home/agmcfarland/flu_project/test/test3/sampleOutputs/IBV_Yamagata_Ref_snpindel'
#opt$samplename <- 'IBV_Victoria_Ref_perfect'
#opt$sampleDir <- '/home/agmcfarland/flu_project/test/test3/sampleOutputs/IBV_Victoria_Ref_perfect'
#troubleshooting inputs sample end

#troubleshooting inputs run start
#opt$softwareDir <- '/home/agmcfarland/flu_project/flu_pipeline'
#opt$report_type <- 'run'
#opt$baseDir <- '/home/agmcfarland/flu_project/test/test4'
#troubleshooting inputs run end

if (opt$report_type=='sample'){
  rmarkdown::render(file.path(opt$softwareDir, 'sample_report.Rmd'),
                    output_file = file.path(opt$sampleDir, paste0(opt$samplename, '.pdf')),
                    params = list(
                      date  = format(Sys.time(), "%Y-%m-%d"),
                      title = paste0('Influenza ', opt$samplename),
                      sampleDir=opt$sampleDir,
                      samplename=opt$samplename))
  
  } else if (opt$report_type=='run') {
      rmarkdown::render(file.path(opt$softwareDir, 'run_report.Rmd'),
                    output_file = file.path(opt$baseDir, 'run_summary.pdf'),
                    params = list(
                      date  = format(Sys.time(), "%Y-%m-%d"),
                      title = paste0('Influenza Run Summary'),
                      baseDir = opt$baseDir))
    } else {
  print('nothing')
}



