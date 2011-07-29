package ru.spbau.bioinf.shift;

import ru.spbau.bioinf.shift.util.ReaderUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

public class Configuration {

    private String proteinDatabaseName;

    private File inputDir;
    private File resultDir;
    private File xmlDir;
    private File xmlSpectrumsDir;
    private File xmlProteinsDir;

    private File datasetDir;

    private ScoringFunction scoringFunction;

    public Configuration(String args[]) {
        String dataset = "data/set8";
        proteinDatabaseName = "prot.fasta";

        if (args != null) {
            if (args.length > 0) {
                dataset = args[0];
            }

            if (args.length > 1) {
                proteinDatabaseName = args[1];
            }
        }
        datasetDir = new File(dataset);
        inputDir = createDir("input");
        resultDir = createDir("result");
        xmlDir = createDir("xml");
        xmlSpectrumsDir = createDir(xmlDir, "spectrums");
        xmlProteinsDir = createDir(xmlDir, "proteins");

        createDir("html");
        scoringFunction = new SharedPeaksScoringFunction();
    }

    public ScoringFunction getScoringFunction() {
        return scoringFunction;
    }

    private File createDir(String name) {
        File dir = new File(datasetDir, name);
        dir.mkdirs();
        return dir;
    }

    private File createDir(File parent, String name) {
        File dir = new File(parent, name);
        dir.mkdirs();
        return dir;
    }

    public File getSpectrumsFile() {
        return new File(inputDir, "input_data");
    }

    public File getProteinDatabaseFile() {
        return new File(inputDir, proteinDatabaseName);
    }

    public File getMatchFile() {
        return new File(resultDir, "match.txt");
    }

    public File getSecondMatchFile() {
        return new File(resultDir, "match2.txt");
    }

    public File getPairsMatchFile() {
        return new File(resultDir, "match_pairs.txt");
    }

    public File getSpectrumXmlFile(Spectrum spectrum) {
        return new File(xmlSpectrumsDir, "spectrum" + spectrum.getId() + ".xml");
    }

    public File getProteinXmlFile(Protein protein) {
        return new File(xmlProteinsDir, "protein" + protein.getProteinId() + ".xml");
    }

    public File getProteinXmlFile() {
        return new File(xmlDir, "proteins.xml");
    }

    public List<Spectrum> getSpectrums() throws IOException {
        BufferedReader input = ReaderUtil.getBufferedReader(getSpectrumsFile());

        Properties properties;
        List<Spectrum> spectrums = new ArrayList<Spectrum>();
        while ((properties = ReaderUtil.readPropertiesUntil(input, "PRECURSOR_MASS")).size() > 0) {
            Spectrum spectrum = new Spectrum(properties, input);
            spectrums.add(spectrum);
        }
        return spectrums;
    }

    public List<Protein> getProteins() throws Exception {
        ProteinDatabaseReader databaseReader = new ProteinDatabaseReader(getProteinDatabaseFile());
        return databaseReader.getProteins();
    }
}
