package ru.spbau.bioinf.shift;

import java.util.Arrays;
import java.util.List;

public enum Modifications {

    BONLY(){
        @Override
        List<Double> getModifications(double p, double precursorMass) {
            return Arrays.asList(p);
        }
    },

    CORE(){
        @Override
        List<Double> getModifications(double p, double precursorMass) {
            double r = precursorMass - p - Y_SHIFT;
            return Arrays.asList(p, r);
        }
    },
    ADDITIONAL(){
        @Override
        List<Double> getModifications(double p, double precursorMass) {
            double r = precursorMass - p - Y_SHIFT;
            return Arrays.asList(p-1, p+1, p + Consts.CO, p + Consts.WATER, r-1, r+1, r - Consts.AMMONIA, r - Consts.WATER, r + Consts.WATER);
        }
    };

    private static final double Y_SHIFT = 0.035;

    abstract List<Double> getModifications(double p, double precursorMass);

}
