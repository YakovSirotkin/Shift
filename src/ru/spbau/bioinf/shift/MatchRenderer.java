package ru.spbau.bioinf.shift;

import org.apache.log4j.Logger;
import org.jdom.Document;
import org.jdom.Element;
import ru.spbau.bioinf.shift.util.ReaderUtil;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;

public class MatchRenderer {

    private static final Logger log = Logger.getLogger(MatchRenderer.class);

    public static void main(String[] args) throws Exception {
        String dataset = "data/set8";
        String proteinDatabaseName = "prot.fasta";

        if (args.length > 0) {
            dataset = args[0];
        }

        if (args.length > 1) {
            proteinDatabaseName = args[1];
        }

        MatchRenderer renderer = new MatchRenderer();


        renderer.render(new File(dataset), proteinDatabaseName);
    }

    public void render(File datasetDir, String proteinDatabaseFilename) throws Exception {
        log.debug("star result processing");
        File msinputDir = new File(datasetDir, "msinput");
        BufferedReader input = ReaderUtil.getBufferedReader(new File(msinputDir, "input_data"));
        List<Spectrum> spectrums = new ArrayList<Spectrum>();
        File outputDir = new File(datasetDir, "xml");
        File msoutputDir = new File(datasetDir, "msoutput");
        File matchFile = new File(msoutputDir, "match.txt");

        File xmlDir = new File(outputDir, "spectrums");
        xmlDir.mkdirs();

        Properties properties;

        while ((properties = ReaderUtil.readPropertiesUntil(input, "PRECURSOR_MASS")).size() > 0) {
            Spectrum spectrum = new Spectrum(properties, input);
            spectrums.add(spectrum);
        }
        log.debug("spectrums data loaded");
        ProteinDatabaseReader fastaReader = new ProteinDatabaseReader(new File(msinputDir, proteinDatabaseFilename));
        List<Protein> proteins = fastaReader.getProteins();

        BufferedReader matchReader = ReaderUtil.getBufferedReader(matchFile);
        HashMap<Integer, List<Integer>> res = new HashMap<Integer, List<Integer>>();
        String s;
        while ((s = matchReader.readLine()) != null) {
            String[] data = s.split(" ");
            int spectrumId = Integer.parseInt(data[0]);
            int proteinId = Integer.parseInt(data[1]);
            if (!res.containsKey(spectrumId)) {
                res.put(spectrumId, new ArrayList<Integer>());
            }
            res.get(spectrumId).add(proteinId);
        }

        for (Map.Entry<Integer, List<Integer>> entry : res.entrySet()) {
            int spectrumId = entry.getKey();
            Spectrum spectrum = spectrums.get(spectrumId);
            List<Integer> proteinIds = entry.getValue();
            Document doc = new Document();
            Element root = new Element("spectrum");
            doc.setRootElement(root);
            spectrum.addToXml(root, new HashMap<Integer, List<Break>>());
            Element matches = new Element("matches");
            for (int proteinId  : proteinIds) {
                Protein protein  = proteins.get(proteinId);
                double[] sd = spectrum.getData();
                double[] pd = protein.getSpectrum();
                Map<String,List<Double>> positions = ProteinFinder.getPositions(spectrum);
                List<Double> shifts = ProteinFinder.getShifts(protein, positions);
                Match m = new Match(protein);
                for (double shift : shifts) {
                    m.addShift(shift, ProteinFinder.getScore(sd, pd, shift));
                }
                matches.addContent(m.toXml());
            }
            root.addContent(matches);

            XmlUtil.saveXml(doc, outputDir + "/spectrums/spectrum" + spectrum.getId() + ".xml");
            log.debug("Match for spectrum " + spectrum.getId() + " saved");

        }
    }

}
