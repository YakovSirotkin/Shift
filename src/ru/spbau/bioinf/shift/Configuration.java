package ru.spbau.bioinf.shift;

import ru.spbau.bioinf.shift.util.ReaderUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;

public class Configuration {

    private File proteinDatabase;

    private File inputDir;
    private File resultDir;
    private File xmlDir;
    private File xmlSpectrumsDir;
    private File xmlProteinsDir;

    private File datasetDir;

    private ScoringFunction scoringFunction;
    private File inputData;

    public Configuration(String args[]) {
        String dataset = "data/set8";
        if (args != null) {
            if (args.length > 0) {
                dataset = args[0];
            }
        }
        datasetDir = new File(dataset);
        inputDir = createDir("input");
        File[] proteinDatabases = inputDir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(".fasta");
            }
        });
        if (proteinDatabases.length == 1) {
            proteinDatabase = proteinDatabases[0];
        } else {
            proteinDatabase = new File(inputDir, args[1]);
        }

        inputData = new File(inputDir, "input_data");

        if (!inputData.exists()) {
            File[] msalignFiles = inputDir.listFiles(new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.endsWith(".msalign");
                }
            });
            if (msalignFiles.length == 1) {
                inputData = msalignFiles[0];
            }
        }


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
        return inputData;
    }

    public File getProteinDatabaseFile() {
        return proteinDatabase;
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

    public File getSignalPeptidesFile() {
        return new File(resultDir, "sp.fasta");
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

    public Map<Integer, Spectrum> getSpectrums() throws IOException {
        File scanDir = new File(inputDir, "env");
        if (scanDir.exists()) {
            return getScans(scanDir);
        }

        BufferedReader input = ReaderUtil.getBufferedReader(getSpectrumsFile());

        Properties properties;
        Map<Integer, Spectrum>  spectrums = new HashMap<Integer, Spectrum>();
        while ((properties = ReaderUtil.readPropertiesUntil(input, "PRECURSOR_MASS")).size() > 0) {
            Spectrum spectrum = new Spectrum(properties, input);
            spectrums.put(spectrum.getId(), spectrum);
        }
        return spectrums;
    }

    public List<MSAlignDiff.MsMatch> getMSAlignResults() throws IOException {
        File resultFile = new File(resultDir, "target_result_detail");

        BufferedReader input = ReaderUtil.getBufferedReader(resultFile);

        Properties properties;
        List<MSAlignDiff.MsMatch> ans = new ArrayList<MSAlignDiff.MsMatch>();
        while ((properties = ReaderUtil.readPropertiesUntil(input, "END PRSM")).size() > 0) {
            ans.add(new MSAlignDiff.MsMatch(properties));
        }
        return ans;
    }

    private Map<Integer, Spectrum> getScans(File scanDir) throws IOException {
        File[] files = scanDir.listFiles(new FileFilter() {
            public boolean accept(File pathname) {
                return pathname.getName().endsWith(".env");
            }
        });

        Map<Integer, Spectrum> spectrums = new HashMap<Integer, Spectrum>();
        for (File file : files) {
            BufferedReader input = ReaderUtil.getBufferedReader(file);

            Properties properties = ReaderUtil.readPropertiesUntil(input, "BEGIN ENVELOPE");
            String fileName = file.getName();
            int id = Integer.parseInt(fileName.substring(fileName.lastIndexOf("_") + 1, fileName.lastIndexOf(".")));
            Spectrum spectrum = new Spectrum(properties, input, id);
            spectrums.put(spectrum.getId(), spectrum);
        }
        return spectrums;
    }


    public List<Protein> getProteins() throws Exception {
        ProteinDatabaseReader databaseReader = new ProteinDatabaseReader(getProteinDatabaseFile());
        return databaseReader.getProteins();
    }
}
