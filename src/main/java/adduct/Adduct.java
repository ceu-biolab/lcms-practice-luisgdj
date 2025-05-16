package adduct;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/* Adduct: product of a direct addition of two or more distinct molecules,
 * resulting in a single reaction product containing all atoms of all components.
 */
public abstract class Adduct {

    /**
     * Calculate the mass to search depending on the adduct hypothesis.
     * Supports single and multiple charges, as well as multimers (e.g., [2M+Na]+).
     *
     * @param mz mz
     * @param adduct adduct notation ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the monoisotopic mass of the experimental mass mz with the adduct @param adduct
     */
    public static Double getMonoisotopicMassFromMZ(Double mz, Map.Entry<String, Double> adduct) {

        Double monoisotopicMass = null;
        Double adductMass = adduct.getValue();
        // !! TODO METHOD
        // !! TODO Create the necessary regex to obtain the multimer (number before the M) and the charge (number before the + or - (if no number, the charge is 1).
        int multimer = extractMultimer(adduct.getKey());
        int charge = extractCharge(adduct.getKey());

        // Case 1: Single charge, no multimer (m/z = M +- adductMass).
        if (charge == 1 && multimer == 1) {
            monoisotopicMass = mz + adductMass;
        }
        // Case 2: Multiple charges, no multimer (mz = M/charge +- adductMass).
        else if (charge > 1 && multimer == 1) {
            monoisotopicMass = (mz + adductMass) * charge;
        }
        // Case 3: Multimer, single charge (mz = M * numberOfMultimer +- adductMass).
        else if (charge == 1 && multimer > 1) {
            monoisotopicMass = (mz + adductMass) / multimer;
        }
        // Case 4: Multimer with multiple charges (mz = (M * numberOfMultimer +- adductMass) / charge).
        else {
            monoisotopicMass = ((mz + adductMass) * charge) / multimer;
        }
        return monoisotopicMass;
    }

    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     * Supports single and multiple charges, as well as multimers (e.g., [2M+Na]+).
     *
     * @param monoisotopicMass The neutral mass of the molecule (M).
     * @param adduct adduct notation (e.g. "[M+H]+", "[2M+Na]+", "[M+2H]2+").
     *
     * @return The expected m/z for the adduct.
     */
    public static Double getMZFromMonoisotopicMass(Double monoisotopicMass, Map.Entry<String, Double> adduct) {

        Double mz = null;
        Double adductMass = adduct.getValue();
        // !! TODO METHOD
        // !! TODO Create the necessary regex to obtain the multimer (number before the M) and the charge (number before the + or - (if no number, the charge is 1).
        int multimer = extractMultimer(adduct.getKey());
        int charge = extractCharge(adduct.getKey());

        // Case 1: Single charge, no multimer (m/z = M +- adductMass).
        if (charge == 1 && multimer == 1) {
            mz = monoisotopicMass - adductMass;
        }
        // Case 2: Multiple charges, no multimer (mz = M/charge +- adductMass).
        else if (charge > 1 && multimer == 1) {
            mz = (monoisotopicMass / charge) - adductMass;
        }
        // Case 3: Multimer, single charge (mz = M * numberOfMultimer +- adductMass).
        else if (charge == 1 && multimer > 1) {
            mz = (monoisotopicMass * multimer) - adductMass;
        }
        // Case 4: Multimer with multiple charges (mz = (M * numberOfMultimer +- adductMass) / charge).
        else {
            mz = ((monoisotopicMass * multimer) - adductMass) / charge;
        }
        return mz;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param theoreticalMass Theoretical mass of the compound
     */
    public static int calculatePPMIncrement(Double experimentalMass, Double theoreticalMass) {
        int ppmIncrement;
        ppmIncrement = (int) Math.round(Math.abs((experimentalMass - theoreticalMass)*1000000 / theoreticalMass));
        return ppmIncrement;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param ppm ppm of tolerance
     */
    public static double calculateDeltaPPM(Double experimentalMass, int ppm) {
        double deltaPPM;
        // deltaPPM = Math.round(Math.abs((experimentalMass * ppm) / 1000000));
        deltaPPM = Math.abs((experimentalMass * ppm) / 1000000);
        return deltaPPM;
    }

    /**
     * Extracts the multimer number from the adduct string using a regex pattern.
     * The multimer refers to the number of molecules in the adduct. If no number is present,
     * the method defaults to 1.
     *
     * @param adduct The adduct string in standard format (e.g. "[2M+H]+").
     * @return The multimer count as an integer (default 1 if not found).
     */
    private static int extractMultimer(String adduct) {

        // Pattern looks for a number before 'M', e.g. [2M or [M
        Matcher m = Pattern.compile("\\[([0-9]*)M").matcher(adduct);
        if (m.find()) {
            String num = m.group(1); // Group 1 captures the digits before 'M'
            return num.isEmpty() ? 1 : Integer.parseInt(num);
        }
        return 1; // Default: no multimer specified
    }

    /**
     * Extracts the charge value from the adduct string using a regex pattern.
     * The charge is found as a number directly before the '+' or '-' at the end of the adduct.
     * If no number is present, the method assumes a default charge of 1.
     *
     * @param adduct The adduct string in standard format (e.g. "[M+2H]2+").
     * @return The charge value (default 1 if not explicitly stated).
     */
    private static int extractCharge(String adduct) {

        // Match the last digit(s) followed by + or − before the closing bracket
        Matcher m = Pattern.compile("([0-9]*)([+-])\\]?$").matcher(adduct);
        if (m.find()) {
            String num = m.group(1);  // May be empty
            return num.isEmpty() ? 1 : Integer.parseInt(num);
        }
        return 1; // Default if no explicit charge found
    }
}