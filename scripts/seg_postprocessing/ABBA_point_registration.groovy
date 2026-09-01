/**
 * ABBA: pixel coordinates -> Allen atlas region name, for a list of points on one sample.
 *
 * REQUIREMENTS
 *   1. The image is open in QuPath.
 *   2. ABBA registration has been exported to the QuPath project.
 *   3. Warped atlas annotations are loaded into the image.
 *
 * If not yet loaded, this script tries to load them automatically.
 */

import static qupath.lib.gui.scripting.QPEx.*
import qupath.ext.biop.abba.AtlasTools

String getAllenRegionId(def ann) {
    // 1) Try metadata keys first
    def md = ann.getMetadata()
    if (md != null) {
        for (key in ["id", "ID", "allen_id", "atlas_id", "structure_id"]) {
            def val = md.get(key)
            if (val != null && val.toString().isInteger()) {
                return val.toString()
            }
        }
    }

    // 2) Fallback: extract a numeric token from the annotation name
    def name = ann.getName()
    if (name != null) {
        def m = (name =~ /(\\d+)/)
        if (m.find()) {
            return m.group(1)
        }
    }

    return ""
}

// ======================== EDIT ME ========================
def IN_CSV  = "path/to/points.csv"
def OUT_CSV = "path/to/points_mapped.tsv" //.tsv due to existance of ',' in Allen Atlas names
def REGION_LABEL_PROPERTY = "name"     // "name" for full Allen names, "acronym" for acronyms
def REGION_ID_PROPERTY = "id"
def SPLIT_HEMISPHERES = true   // usually true
def EXCLUDE_LEFT_RIGHT = true  // skip hemisphere container annotations as final region labels
// =========================================================

def imageData = getCurrentImageData()

println "image: " + imageData.getServer().getMetadata().getName()

def loadAtlasAnnotationsWithProperty = { imgData, propertyName, splitHemispheres, excludeLeftRight ->
    removeObjects(getAnnotationObjects(), false)
    AtlasTools.loadWarpedAtlasAnnotations(imgData, propertyName, splitHemispheres)

    def anns = getAnnotationObjects().findAll { a ->
        def n = a.getName()
        if (n == null || n.trim().isEmpty()) return false
        if (excludeLeftRight && (n == "Left" || n == "Right")) return false
        return a.getROI() != null
    }.sort { a, b -> a.getROI().getArea() <=> b.getROI().getArea() }

    return anns
}
def annsLabel = loadAtlasAnnotationsWithProperty(
    imageData,
    REGION_LABEL_PROPERTY,
    SPLIT_HEMISPHERES,
    EXCLUDE_LEFT_RIGHT
)
def annsId = loadAtlasAnnotationsWithProperty(
    imageData,
    REGION_ID_PROPERTY,
    SPLIT_HEMISPHERES,
    EXCLUDE_LEFT_RIGHT
)

println "Label annotations: ${annsLabel.size()}"
println "ID annotations: ${annsId.size()}"

def inFile = new File(IN_CSV)
if (!inFile.exists()) {
    println "ERROR: input file not found: " + IN_CSV
    return
}

def nPoints = 0
inFile.eachLine { line, i ->
    if (i == 1) return   // skip header
    if (line != null && !line.trim().isEmpty()) nPoints++
}
println "Number of points to map: ${nPoints}"

int n = 0
int matched = 0
int unmatched = 0

new File(OUT_CSV).withWriter { writer ->
    inFile.withReader { reader ->
        def firstLine = reader.readLine()
        if (firstLine == null) {
            println "ERROR: empty input CSV"
            return
        }

        def header = firstLine.split(",").toList()
        int iId = header.indexOf("cell_id")
        int iPx = header.indexOf("x")
        int iPy = header.indexOf("y")

        if (iId < 0 || iPx < 0 || iPy < 0) {
            println "ERROR: CSV must have cell_id, px, py columns. Found: " + header
            return
        }

        writer.writeLine((header + ["atlas_region", "atlas_region_id"]).join("\t"))

        String line
        while ((line = reader.readLine()) != null) {
            if (line.trim().isEmpty()) continue

            def f = line.split(",")
            def id = f[iId]
            double px = Double.parseDouble(f[iPx])
            double py = Double.parseDouble(f[iPy])

            def hitsLabel = annsLabel.findAll { ann -> ann.getROI().contains(px, py) }
            def hitsId    = annsId.findAll    { ann -> ann.getROI().contains(px, py) }

            String atlasRegion = ""
            String atlasRegionId = ""

            if (!hitsLabel.isEmpty()) {
                atlasRegion = hitsLabel[0].getName() ?: ""
            }

            if (!hitsId.isEmpty()) {
                atlasRegionId = hitsId[0].getName() ?: ""
            }

            if (atlasRegion || atlasRegionId) {
                matched++
            } else {
                unmatched++
            }

            writer.writeLine((f + [atlasRegion, atlasRegionId]).join("\t"))

            n++
            if (n % 10000 == 0) {
                println "  ${n} points processed..."
            }
        }
    }
}

println "wrote ${n} rows -> ${OUT_CSV}"
println "matched: ${matched}"
println "unmatched: ${unmatched}"
