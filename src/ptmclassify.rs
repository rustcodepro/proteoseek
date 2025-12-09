use smartcore::linalg::basic::matrix::DenseMatrix;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
/*
Gaurav Sablok
codeprog@icloud.com
*/

type VECALLOCATE = (String, usize);
type VECCOMMONALLOCATION = (String, usize, usize);
type LASTCOMMONALLOCATION = (String, usize, usize, usize);
type DENSEMATRIX = (DenseMatrix<f32>, DenseMatrix<f32>, Vec<i32>);

pub fn classify(
    pathfile: &str,
    pathannotate: &str,
    pathkmer: &str,
    threasholdvalue: &str,
    predictionfastafile: &str,
    predictionpeak: &str,
) -> Result<DENSEMATRIX, Box<dyn Error>> {
    let mut peakdata: Vec<i32> = Vec::new();
    let peakopen = File::open(pathannotate).expect("file not present");
    let peakread = BufReader::new(peakopen);
    for i in peakread.lines() {
        let line = i.expect("file not present");
        let value = line.parse::<i32>().unwrap();
        peakdata.push(value);
    }
    let mut pathmervec: Vec<String> = Vec::new();
    let pathmervecopen = File::open(pathkmer).expect("file not present");
    let pathmervecread = BufReader::new(pathmervecopen);
    for i in pathmervecread.lines() {
        let line = i.expect("file not present");
        pathmervec.push(line);
    }
    let mut pathfastavec: Vec<String> = Vec::new();
    let pathfastaopen = File::open(pathfile).expect("file not present");
    let pathfastaread = BufReader::new(pathfastaopen);
    for i in pathfastaread.lines() {
        let line = i.expect("line not present");
        pathfastavec.push(line);
    }
    let mut pathkmer: Vec<Vec<String>> = Vec::new();
    for i in pathfastavec.iter() {
        let insertkmer = i
            .as_bytes()
            .windows(3)
            .map(|x| std::str::from_utf8(x).unwrap())
            .map(|x| x.to_string())
            .collect::<Vec<_>>();
        pathkmer.push(insertkmer);
    }
    let mut ptmvec: Vec<VECALLOCATE> = Vec::new();
    for i in pathkmer.iter() {
        for val in i.iter() {
            for kmervec in pathmervec.iter() {
                let mut count = 0usize;
                if val.clone() == kmervec.clone() {
                    count += 1usize
                }
                ptmvec.push((val.clone(), count));
            }
        }
    }

    let mut commonallocation: Vec<VECCOMMONALLOCATION> = Vec::new();
    for i in ptmvec.iter() {
        for val in peakdata.iter() {
            let valueinsert: (String, usize, usize) = (i.0.clone(), i.1, *val as usize);
            commonallocation.push(valueinsert);
        }
    }

    /*
     association criteria: algorithmic view: Number of the kmers occurences from the previous
     sequences using the PTM generate and then the peak data clustering.
     if the value is there and more than the threshold
     bumping up the cluster and if the value is less then lowering the cluster, so that the
     binary classification can be configured.
    */

    let mut commonallocation_joined_array: Vec<LASTCOMMONALLOCATION> = Vec::new();
    let mut classification_label: Vec<i32> = Vec::new();
    for i in commonallocation.iter() {
        if i.2 <= threasholdvalue.parse::<usize>().unwrap() {
            let value: (String, usize, usize, usize) = (i.0.clone(), i.1, i.2, (i.1 + i.2));
            commonallocation_joined_array.push(value);
            classification_label.push(0);
        }
        if i.2 > threasholdvalue.parse::<usize>().unwrap() {
            let value: (String, usize, usize, usize) = (i.0.clone(), i.1, i.2, (i.1 - i.2));
            commonallocation_joined_array.push(value);
            classification_label.push(1);
        }
    }

    /*
    Making a uniform padded sequence tensor for the classification:
    */

    let sequence = commonallocation
        .iter()
        .map(|x| x.0.clone())
        .collect::<Vec<String>>();
    let sequencematrix = sequencevector(sequence).unwrap();
    let mut finalsequencematrix: Vec<Vec<f32>> = Vec::new();
    for i in sequencematrix.iter() {
        for val in commonallocation.iter() {
            let value = i;
            let mut valueappend: Vec<f32> = vec![val.1 as f32, val.2 as f32, 0 as f32, 0 as f32];
            value.clone().append(&mut valueappend);
            finalsequencematrix.push(value.clone());
        }
    }

    let predictionvalue: (Vec<String>, Vec<f32>) =
        prediction(predictionfastafile, predictionpeak).unwrap();
    let mut finalpredictionvec: Vec<Vec<f32>> = Vec::new();
    let predictionpad = sequencevector(predictionvalue.0).unwrap();
    for i in predictionpad.iter() {
        for val in predictionvalue.1.iter() {
            let value = i;
            let valueappend: Vec<f32> = vec![*val, 0 as f32, 0 as f32, 0 as f32];
            value.clone().append(&mut valueappend.clone());
            finalpredictionvec.push(value.clone());
        }
    }

    let matrixreturn: DenseMatrix<f32> = DenseMatrix::from_2d_vec(&finalsequencematrix).unwrap();
    let predictionmatrix: DenseMatrix<f32> = DenseMatrix::from_2d_vec(&finalpredictionvec).unwrap();

    let densematrix: DENSEMATRIX = (matrixreturn, predictionmatrix, classification_label);

    Ok(densematrix)
}

pub fn sequencevector(input: Vec<String>) -> Result<Vec<Vec<f32>>, Box<dyn Error>> {
    let valueinsert = input.clone();
    let mut vecfinal: Vec<Vec<f32>> = Vec::new();
    for i in valueinsert.iter() {
        let stringclone = i.chars().collect::<Vec<char>>();
        let mut vecinsert: Vec<Vec<f32>> = Vec::new();
        for i in stringclone.iter() {
            match i {
                'A' => vecinsert.push(vec![1.0, 0.0, 0.0, 0.0]),
                'T' => vecinsert.push(vec![0.0, 1.0, 0.0, 0.0]),
                'G' => vecinsert.push(vec![0.0, 0.0, 1.0, 0.0]),
                'C' => vecinsert.push(vec![0.0, 0.0, 0.0, 1.0]),
                'N' => vecinsert.push(vec![1.0, 1.0, 1.0, 1.0]),
                _ => continue,
            }
        }
        let unpackedvec = vecinsert.iter().flatten().cloned().collect::<Vec<f32>>();
        vecfinal.push(unpackedvec);
    }
    let maxlength: usize = *vecfinal
        .iter()
        .map(|x| x.len())
        .collect::<Vec<_>>()
        .iter()
        .max()
        .unwrap();
    let mut finalpaddedvec: Vec<Vec<f32>> = Vec::new();
    for i in vecfinal.iter() {
        let mut vectinsert: Vec<Vec<f32>> = Vec::new();
        let ilength = i.len();
        if ilength == maxlength {
            finalpaddedvec.push(i.clone());
        }
        if ilength != maxlength {
            let value = maxlength - ilength;
            for _ in 0..value {
                vectinsert.push(vec![0.0, 0.0, 0.0, 0.0]);
            }
        }
        let mut unpackedvecinsert = vectinsert.iter().flatten().cloned().collect::<Vec<f32>>();
        let mut initialvec: Vec<f32> = i.clone();
        initialvec.append(&mut unpackedvecinsert);
        finalpaddedvec.push(initialvec);
    }
    Ok(finalpaddedvec)
}

pub fn prediction(
    pathfile: &str,
    peakdata: &str,
) -> Result<(Vec<String>, Vec<f32>), Box<dyn Error>> {
    let fileopen = File::open(pathfile).expect("file not present");
    let fileread = BufReader::new(fileopen);
    let mut sequencevec: Vec<String> = Vec::new();
    for i in fileread.lines() {
        let line = i.expect("line not present");
        if !line.starts_with("#") {
            sequencevec.push(line)
        }
    }
    let mut peakvecdense: Vec<f32> = Vec::new();
    let peakline = File::open(peakdata).expect("line not present");
    let peakread = BufReader::new(peakline);
    for i in peakread.lines() {
        let line = i.expect("line not present");
        peakvecdense.push(line.parse::<f32>().unwrap());
    }
    Ok((sequencevec, peakvecdense))
}
