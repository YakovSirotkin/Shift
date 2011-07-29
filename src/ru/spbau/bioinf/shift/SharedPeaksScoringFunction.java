package ru.spbau.bioinf.shift;

public class SharedPeaksScoringFunction implements ScoringFunction {

    public double getScore(Spectrum spectrum, Protein protein, double shift) {
        double[] sd = spectrum.getData();
        double[] pd = protein.getSpectrum();

        int score = 0;
        int i = 0;
        int j = 0;
        do {
            double diff = sd[j] - pd[i] + shift;
            if (diff < -0.1) {
                j++;
            } else if (diff > 0.1) {
                i++;
            } else {
                score++;
                i++;
                j++;
            }
        } while (i < pd.length && j < sd.length);
        return score;
    }
}
