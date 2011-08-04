package ru.spbau.bioinf.shift;

public class ExtendedSharedPeaksScoringFunction extends SharedPeaksScoringFunction {
    @Override
    public double getScore(Spectrum spectrum, Protein protein, double shift) {
        double a = super.getScore(spectrum, protein, shift);
        double b = getScore(spectrum.getAdditionalSpectrum(), protein, shift);
        return a + b/10;
    }
}
