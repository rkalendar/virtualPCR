import java.io.IOException;
import java.nio.file.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class OutputPath {

    public enum Type { FILE, DIRECTORY }
    public record Target(Path path, Type type) {}

    public static Target parse(String raw) {
        String s = stripQuotes(raw.trim());
        s = expandTilde(s);
        s = expandEnvVars(s);

        Path p = Paths.get(s).normalize().toAbsolutePath();

        // If it exists, decide by actual type:
        if (Files.exists(p)) {
            return new Target(p, Files.isDirectory(p) ? Type.DIRECTORY : Type.FILE);
        }

        // Heuristics for non-existent paths:
        boolean endsWithSep = s.endsWith("/") || s.endsWith("\\");
        boolean looksLikeFileName =
                p.getFileName() != null &&
                p.getFileName().toString().matches(".*\\.[^.]+$"); // has an extension like ".txt", ".fa", etc.

        Type t = endsWithSep ? Type.DIRECTORY : (looksLikeFileName ? Type.FILE : Type.DIRECTORY);
        return new Target(p, t);
    }

    public static void ensureParentExists(Target t) throws IOException {
        if (t.type == Type.DIRECTORY) {
            Files.createDirectories(t.path);
        } else {
            Path parent = t.path.getParent();
            if (parent != null) Files.createDirectories(parent);
        }
    }

    private static String stripQuotes(String s) {
        if (s.length() >= 2 && ((s.startsWith("\"") && s.endsWith("\"")) || (s.startsWith("'") && s.endsWith("'")))) {
            return s.substring(1, s.length() - 1);
        }
        return s;
    }

    private static String expandTilde(String s) {
        if (s.startsWith("~/") || s.startsWith("~\\")) {
            return System.getProperty("user.home") + s.substring(1);
        }
        return s;
        // (You can handle bare "~" if desired.)
    }

    private static String expandEnvVars(String s) {
        // ${VAR} style
        Pattern braced = Pattern.compile("\\$\\{([^}]+)}");
        Matcher m = braced.matcher(s);
        StringBuffer sb = new StringBuffer();
        while (m.find()) {
            String val = System.getenv().getOrDefault(m.group(1), "");
            m.appendReplacement(sb, Matcher.quoteReplacement(val));
        }
        m.appendTail(sb);
        s = sb.toString();

        // %VAR% style (Windows)
        Pattern percent = Pattern.compile("%([^%]+)%");
        m = percent.matcher(s);
        sb = new StringBuffer();
        while (m.find()) {
            String val = System.getenv().getOrDefault(m.group(1), "");
            m.appendReplacement(sb, Matcher.quoteReplacement(val));
        }
        m.appendTail(sb);
        return sb.toString();
    }
}