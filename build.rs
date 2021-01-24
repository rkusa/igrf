use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

struct FieldModel {
    name: &'static str,
    /// Epoch of model.
    epoch: f64,
    /// Min year of model.
    yrmin: f64,
    /// Max year of model.
    yrmax: f64,
    /// Minimum height of each model.
    altmin: f64,
    /// Maximum height of each model.
    altmax: f64,
    /// Main field coefficient.
    max1: usize,
    /// Secular variation coefficient.
    max2: usize,
    /// Acceleration coefficient.
    max3: usize,
    /// Schmidt quasi-normal internal spherical harmonic coeff.
    gh1: Vec<f64>,
    /// Schmidt quasi-normal internal spherical harmonic coeff.
    gh2: Vec<f64>,
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    let out_dir = env::var_os("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("data.rs");

    let mut out = BufWriter::new(File::create(&dest_path).unwrap());

    let data = include_str!("./data/IGRF13.COF");
    let mut lines = data.lines().peekable();

    let mut models = Vec::new();

    while let Some(line) = lines.next() {
        let mut parts = line.split_ascii_whitespace();
        let mut model = FieldModel {
            name: parts.next().expect("name part"),
            epoch: FromStr::from_str(parts.next().expect("epoch part")).unwrap(),
            max1: FromStr::from_str(parts.next().expect("epoch part")).unwrap(),
            max2: FromStr::from_str(parts.next().expect("epoch part")).unwrap(),
            max3: FromStr::from_str(parts.next().expect("epoch part")).unwrap(),
            yrmin: FromStr::from_str(parts.next().expect("yrmin part")).unwrap(),
            yrmax: FromStr::from_str(parts.next().expect("yrmax part")).unwrap(),
            altmin: FromStr::from_str(parts.next().expect("altmin part")).unwrap(),
            altmax: FromStr::from_str(parts.next().expect("altmax part")).unwrap(),
            gh1: Vec::new(),
            gh2: Vec::new(),
        };

        let mut gh1 = Vec::with_capacity(model.max1);
        let mut gh2 = Vec::with_capacity(model.max2);

        gh1.push(0.0);
        gh2.push(0.0);

        'outer: loop {
            if let Some(line) = lines.peek() {
                if line.starts_with("   ") {
                    break 'outer;
                }
            }

            if let Some(line) = lines.next() {
                let mut parts = line.split_ascii_whitespace();
                let n = usize::from_str(parts.next().expect("n part")).unwrap();
                let m = usize::from_str(parts.next().expect("m part")).unwrap();
                let g = f64::from_str(parts.next().expect("g part")).unwrap();
                let hh = f64::from_str(parts.next().expect("h part")).unwrap();
                let g2 = f64::from_str(parts.next().expect("g2 part")).unwrap();
                let hh2 = f64::from_str(parts.next().expect("h2 part")).unwrap();

                gh1.push(g);
                if m != 0 {
                    gh1.push(hh);
                }

                if n <= model.max2 {
                    gh2.push(g2);
                    if m != 0 {
                        gh2.push(hh2);
                    }
                }
            } else {
                break;
            }
        }

        model.gh1 = gh1;
        model.gh2 = gh2;
        models.push(model);
    }

    for i in 0..models.len() {
        if models[i].max2 == 0 {
            models[i].gh2 = models[i + 1].gh1.clone();
        }
    }

    write!(out, "const MODELS: &[FieldModel] = &[").unwrap();

    for model in models {
        write!(
            out,
            r#"
                FieldModel {{
                    name: "{}",
                    epoch: {}f64,
                    yrmin: {}f64,
                    yrmax: {}f64,
                    altmin: {}f64,
                    altmax: {}f64,
                    max1: {},
                    max2: {},
                    max3: {},
                    gh1: &[
            "#,
            model.name,
            model.epoch,
            model.yrmin,
            model.yrmax,
            model.altmin,
            model.altmax,
            model.max1,
            model.max2,
            model.max3,
        )
        .unwrap();

        for v in model.gh1 {
            write!(out, "{}f64,", v).unwrap();
        }

        write!(
            out,
            r#"
                    ],
                    gh2: &[
            "#,
        )
        .unwrap();

        for v in model.gh2 {
            write!(out, "{}f64,", v).unwrap();
        }

        write!(
            out,
            r#"
                    ],
                }},
            "#,
        )
        .unwrap();
    }

    write!(out, "];").unwrap();
}
