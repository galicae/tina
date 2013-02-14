/******************************************************************************
 * util.CommandLineParser                                                    *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package util;

/**
 * CommandLineParser is a class that you can give your command line
 * arguments and it will parse them for you in an easy way.
 * It isn't thought to be very effective, but could in time easily be
 * changed.
 * @author huberste
 * @date November 11, 2012
 * @version 1.0
 */
public class CommandLineParser {
	
	private String[] arguments;
	
	/**
	 * Constructs a CommandLineParser.
	 * @param args the CommandLine Arguments to be parsed, normally 'args'.
	 */
	public CommandLineParser(String[] args) {
		arguments = new String[args.length];
		for (int i = 0; i < args.length; i++) {
			arguments[i]=args[i];
		}
	}
	
	/**
	 * Checks if specific true/false argument was given
	 * @param arg the name of the argument
	 * @return true if given, false else
	 */
	public boolean getBoolArg(String arg) {
		for (int i = 0; i < arguments.length; i++) {
			if (arguments[i].equals(arg)) {
				return true;
			}
		}
		return false;
	}
	/**
	 * Checks if specific true/false argument was given
	 * @param arg the name of the argument
	 * @param defaultValue the default value that will be returned if no
	 * valid argument was given
	 * @return true if given, false else
	 */
	public boolean getBoolArg(String arg, boolean defaultValue) {
		boolean result = false;
		for (int i = 0; i < arguments.length; i++) {
			if (arguments[i].equals(arg)) {
				result = true;
				break;
			}
		}
		return result;
	}

	/**
	 * Checks for specific int argument given
	 * @param arg the name of the argument given
	 * @return given integer, 0 if none was given
	 */
	public int getIntArg(String arg) {
		for (int i = 0; i < arguments.length; i++) {
			if (arguments[i].equals(arg)) {
				return Integer.parseInt(arguments[i+1]);
			}
		}
		return 0;
	}
	
	/**
	 * Checks for specific int argument given
	 * @param arg the name of the argument given
	 * @param defaultValue the default value that will be returned if no
	 * valid argument was given
	 * @return given integer, 0 if none was given
	 */
	public int getIntArg(String arg, int defaultValue) {
		int result = defaultValue;
		for (int i = 0; i < arguments.length; i++) {
			if (arguments[i].equals(arg)) {
				result = Integer.parseInt(arguments[i+1]);
				break;
			}
		}
		return result;
	}

	/**
	 * Checks for specific Double argument was given
	 * @param arg the name of the argument given
	 * @return given double, 0.0 if none was given
	 */
	public double getDoubleArg(String arg) {
		for (int i = 0; i < arguments.length; i++) {
			if (arguments[i].equals(arg)) {
				return Double.parseDouble(arguments [i+1]);
			}
		}
		return 0.0;
	}
	
	/**
	 * Checks for specific Double argument was given
	 * @param arg the name of the argument given
	 * @param defaultValue the default value that will be returned if no
	 * valid argument was given
	 * @return given double, 0.0 if none was given
	 */
	public double getDoubleArg(String arg, double defaultValue) {
		double result = defaultValue;
		for (int i = 0; i < arguments.length; i++) {
			if (arguments[i].equals(arg)) {
				result = Double.parseDouble(arguments [i+1]);
				break;
			}
		}
		return result;
	}

	/**
	 * Checks for specific String argument was given
	 * @param arg the name of the argument given
	 * @return given String, null if none was given
	 */
	public String getStringArg(String arg) {
		for (int i = 0; i < arguments.length; i++) {
			if (arguments[i].equals(arg)) {
				return arguments[i+1];
			}
		}
		return null;
	}
	
	/**
	 * Checks for specific String argument was given
	 * @param arg the name of the argument given
	 * @param defaultValue the default value that will be returned if no
	 * valid argument was given
	 * @return given String, default if none was given
	 */
	public String getStringArg(String arg, String defaultValue) {
		String result = defaultValue;
		for (int i = 0; i < arguments.length; i++) {
			if (arguments[i].equals(arg)) {
				result = arguments[i+1];
				break;
			}
		}
		return result;
	}
	
	/**
	 * Checks for specific char argument was given
	 * @param arg the name of the argument given
	 * @return given char, ' ' if none was given
	 */
	public char getCharArg(String arg) {
		for (int i = 0; i < arguments.length; i++) {
			if (arguments[i].equals(arg)) {
				return arguments[i+1].charAt(0);
			}
		}
		return ' ';
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/