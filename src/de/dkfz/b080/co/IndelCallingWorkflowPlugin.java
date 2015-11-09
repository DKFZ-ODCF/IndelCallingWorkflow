package de.dkfz.b080.co;

import de.dkfz.roddy.plugins.BasePlugin;

/**
 */
public class IndelCallingWorkflowPlugin extends BasePlugin {

    public static final String CURRENT_VERSION_STRING = "1.0.161";
    public static final String CURRENT_VERSION_BUILD_DATE = "Wed May 27 12:15:48 CEST 2015";

    @Override
    public String getVersionInfo() {
        return "Roddy plugin: " + this.getClass().getName() + ", V " + CURRENT_VERSION_STRING + " built at " + CURRENT_VERSION_BUILD_DATE;
    }
}
