use clap::{Parser, Subcommand};
#[derive(Debug, Parser)]
#[command(
    name = "proteoseek",
    version = "1.0",
    about = "proteoseek for protein machine learning
       ************************************************
       Gaurav Sablok
       codeprog@icloud.com
      ************************************************"
)]
pub struct CommandParse {
    /// subcommands for the specific actions
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Generate PTM classification
    GeneratePTM {
        /// path to the previously annotated fasta
        fastafile: String,
        /// threads for the analysis
        thread: String,
    },
    /// PTM machine learning
    PTMClassify {
        /// path to the fasta file
        fastafile: String,
        /// path to the previous kmer
        kmerpath: Option<String>,
        /// annotated peak data
        peakdata: String,
        /// threashold for the analysis
        threshold: String,
        /// prediction fasta file
        predfasta: String,
        /// prediction peak
        predpeak: String,
        /// threads for the analysis
        threads: String,
    },
}
