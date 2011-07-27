package ru.spbau.bioinf.shift;

import org.apache.log4j.Logger;
import org.jdom.Document;
import org.jdom.Element;
import ru.spbau.bioinf.shift.util.ReaderUtil;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MatchRenderer {

    private static final Logger log = Logger.getLogger(MatchRenderer.class);

    public static void main(String[] args) throws Exception {
        Configuration config = new Configuration(args);
        MatchRenderer renderer = new MatchRenderer();
        renderer.render(config);
    }

    public void render(Configuration config) throws Exception {
        log.debug("star result processing");
        List<Spectrum> spectrums = config.getSpectrums();

        log.debug("spectrums data loaded");
        List<Protein> proteins = config.getProteins();

        BufferedReader matchReader = ReaderUtil.getBufferedReader(config.getMatchFile());
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
            Element root = spectrum.toXml( new HashMap<Integer, List<Break>>());
            doc.setRootElement(root);
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

            XmlUtil.saveXml(doc, config.getSpectrumXmlFile(spectrum));
            log.debug("Match for spectrum " + spectrum.getId() + " saved");
        }
    }

}
