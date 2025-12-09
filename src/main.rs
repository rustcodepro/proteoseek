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
use smartcore::tree::decision_tree_classifier::DecisionTreeClassifierParameters;
/*
Gaurav Sablok
codeprog@icloud.com
*/

fn main() {
    let standard_font = FIGfont::standard().unwrap();
    let figure = standard_font.convert("proteoseek");
    assert!(figure.is_some());
    println!("{}", figure.unwrap());
    let argparse = CommandParse::parse();
    match &argparse.command {
        Commands::GeneratePTM { fastafile, thread } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let value = ptmgenerate(fastafile).unwrap();
            })
        }
        Commands::PTMClassify {
            fastafile,
            kmerpath,
            peakdata,
            threshold,
            predfasta,
            predpeak,
            threads,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let valudense: DenseMatrix<f32> = classify(
                    fastafile, peakdata, kmerpath, predfasta, predpeak, threshold,
                )
                .unwrap();
                let machinelearning = LogisticRegression::fit(&valudense.0, &valudense.3).unwrap();
                let machinelearningmodel = machinelearning.predict(&valudense.2).unwrap();
                let accuracy_defined = accuracy(&valudense.2, &machinelearningmodel).abs();
                print!("The accuracy of the model is {:?}", accuracy_defined);
            })
        }
    }
}
