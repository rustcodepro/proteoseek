mod args;
mod ptmclassify;
mod ptmgen;
use crate::args::CommandParse;
use crate::args::Commands;
use crate::ptmclassify::classify;
use crate::ptmgen::ptmgenerate;
use clap::Parser;
use smartcore::linalg::basic::matrix::DenseMatrix;
use smartcore::linear::logistic_regression::LogisticRegression;
use smartcore::metrics::accuracy;

/*
Gaurav Sablok
codeprog@icloud.com
*/

fn main() {
    let argparse = CommandParse::parse();
    match &argparse.command {
        Commands::GeneratePTM { fastafile, thread } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let value = ptmgenerate(fastafile).unwrap();
                print!("The command has finished:{:?}", value);
            })
        }
        Commands::PTMClassify {
            fastafile,
            kmerpath,
            peakdata,
            threshold,
            predfasta,
            predpeak,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let valudense: (DenseMatrix<f32>, DenseMatrix<f32>, Vec<i32>) = classify(
                    fastafile, peakdata, kmerpath, predfasta, predpeak, threshold,
                )
                .unwrap();
                let machinelearning =
                    LogisticRegression::fit(&valudense.0, &valudense.2, Default::default())
                        .unwrap();
                let machinelearningmodel = machinelearning.predict(&valudense.1).unwrap();
                let accuracy_defined = accuracy(&valudense.2, &machinelearningmodel).abs();
                print!("The accuracy of the model is {:?}", accuracy_defined);
            })
        }
    }
}
