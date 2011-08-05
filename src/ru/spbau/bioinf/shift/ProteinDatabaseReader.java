package ru.spbau.bioinf.shift;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class ProteinDatabaseReader {

    private ArrayList<Protein> proteins = new ArrayList<Protein>();

    public ArrayList<Protein> getProteins() {
        return proteins;
    }

    public ProteinDatabaseReader(File proteinDatabase)
            throws Exception {

        BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(
                proteinDatabase), "UTF-8"));

        String s;
        String cur = "";
        int proteinId = 0;
        String name = null;
        while ((s = input.readLine()) != null) {
            if (s.startsWith(">")) {
                if (name != null) {
                    proteins.add(new Protein(proteinId, cur, name));
                    proteinId++;
                    cur = "";
                }
                name = s.substring(1);
            } else {
                cur += s.trim();
            }
        }
        proteins.add(new Protein(proteinId, cur, name));
    }
}
