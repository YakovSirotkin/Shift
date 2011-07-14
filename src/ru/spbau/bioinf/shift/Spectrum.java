package ru.spbau.bioinf.shift;

import org.jdom.Element;
import ru.spbau.bioinf.shift.util.ReaderUtil;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Properties;

public class Spectrum {

    private int id;
    private List<Peak> peaks = new ArrayList<Peak>();
    private String scans;
    private double precursorMz;
    private int precursorCharge;
    private double precursorMass;

    private double[] data = null;

    public Spectrum(Properties prop, BufferedReader input) throws IOException {
        id = ReaderUtil.getIntValue(prop, "ID");
        scans = ReaderUtil.getValue(prop, "SCANS");
        precursorMz = ReaderUtil.getDoubleValue(prop, "PRECURSOR_MZ");
        precursorCharge = ReaderUtil.getIntValue(prop, "PRECURSOR_CHARGE");
        precursorMass = ReaderUtil.getDoubleValue(prop, "PRECURSOR_MASS");
        List<String[]> datas = ReaderUtil.readDataUntil(input, "END IONS");
        for (String[] data : datas) {
            peaks.add(new Peak(Float.parseFloat(data[0]), Float.parseFloat(data[1]), Integer.parseInt(data[2])));
        }
    }

    public Spectrum(List<Peak> peaks, double precursorMass) {
        this.precursorMass = precursorMass;
        this.peaks = peaks;
    }

    public int getId() {
        return id;
    }

    public double getPrecursorMass() {
        return precursorMass;
    }

    public List<Peak> getPeaks() {
        return peaks;
    }

    public int getPeakCharge(int peakId) {
        return peaks.get(peakId).getCharge();
    }

    public double[] getData() {
        if (data == null) {
            List<Double> d = new ArrayList<Double>();
            for (Peak peak : peaks) {
                double p = peak.getMonoisotopicMass();
                d.addAll(getModifications(p));
            }
            Collections.sort(d);
            d = ProteinFinder.merge(d);
            data = new double[d.size()];
            for (int i = 0; i < data.length; i++) {
                data[i] = d.get(i);
            }
        }

        return data;
    }

    private List<Double> getModifications(double p) {
        double r = precursorMass - p;
        return Arrays.asList(
                p - 1, p, p + 1, p + Consts.AMMONIA, p + Consts.WATER,
                r - 1, r, r + 1, r - Consts.AMMONIA, r - Consts.WATER);
    }

    public Spectrum getLeft(double[] p, double shift) {
        List<Peak> left = new ArrayList<Peak>();
        for (Peak peak : peaks) {
            List<Double> mods = getModifications(peak.getMonoisotopicMass());
            boolean found = false;
            for (double mod : mods) {
                for (double v : p) {
                    if (Math.abs(mod - v + shift) < 0.1) {
                        found = true;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
            if (!found) {
                left.add(peak);
            }
        }
        return new Spectrum(left, precursorMass);
    }

    public void clearData() {
        data = null;
    }

    public void addToXml(Element parent, Map<Integer, List<Break>> breaks) {
        Element peaksTag = new Element("peaks");
        parent.addContent(peaksTag);
        XmlUtil.addElement(parent, "scans", scans);
        XmlUtil.addElement(parent, "spectrum-id", id);
        XmlUtil.addElement(parent, "precursor-mz", precursorMz);
        XmlUtil.addElement(parent, "precursor-charge", precursorCharge);
        XmlUtil.addElement(parent, "precursor-mass", precursorMass);
        for (int i = 0; i < peaks.size(); i++) {
            Peak peak = peaks.get(i);
            Element peakTag = new Element("peak");
            peaksTag.addContent(peakTag);
            XmlUtil.addElement(peakTag, "peak-id", i);
            peak.addToXml(peakTag, breaks.get(i));
        }
    }
}
