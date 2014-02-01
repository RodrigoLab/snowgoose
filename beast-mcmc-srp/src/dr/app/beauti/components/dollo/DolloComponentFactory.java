package dr.app.beauti.components.dollo;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ComponentOptions;

/**
 * @author Marc Suchard
 * @version $Id: DolloComponentFactory.java 5861 2013-10-03 13:38:22Z rambaut $
 */

public class DolloComponentFactory implements ComponentFactory {

	public static ComponentFactory INSTANCE = new DolloComponentFactory();

	public ComponentGenerator createGenerator(BeautiOptions beautiOptions) {
        return new DolloComponentGenerator(beautiOptions);
	}

	public ComponentOptions createOptions(BeautiOptions beautiOptions) {
        return new DolloComponentOptions(beautiOptions);
	}

}