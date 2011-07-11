package ru.spbau.bioinf.shift;

import java.io.File;

public class Configuration {

    private String proteinDatabaseName;

    private File inputDir;
    private File resultDir;
    private File xmlDir;
    private File xmlSpectrumsDir;


    private File datasetDir;

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

        createDir("html");
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

    public File getSpectrumXmlFile(Spectrum spectrum) {
        return new File(xmlSpectrumsDir, "spectrum" + spectrum.getId() + ".xml");
    }
}
