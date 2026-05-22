FIRFILT EEGLAB plugin
-------------
The FIRfilt EEGLAB plugin is a tool used within the EEGLAB environment. FIRfilt is specifically designed for filtering EEG data using Finite Impulse Response (FIR) filters. Here are the key features and functionalities of the FIRfilt plugin:

* FIR Filtering: FIRfilt provides a straightforward interface for applying FIR filters to EEG data. FIR filters are commonly used due to their stability and linear phase properties.
* Filter Types: Users can create various types of FIR filters, including low-pass, high-pass, band-pass, and band-stop filters. This flexibility allows users to isolate specific frequency bands of interest.
* Design Methods: FIRfilt offers several methods for designing FIR filters, such as the window method, least-squares method, and equiripple method. Each method has its own advantages depending on the specific filtering requirements.
* Graphical Interface: The plugin integrates with the EEGLAB GUI, making it accessible for users who prefer graphical user interfaces for their data processing tasks.
* Command Line Support: For more advanced users, FIRfilt also supports command-line operations, allowing for script-based automation and integration into larger data processing pipelines.
 
See [this page](https://eeglab.org/others/Firfilt_FAQ.html) or the [paper](https://home.uni-leipzig.de/biocog/eprints/widmann_a2015jneuroscimeth250_34.pdf) for additional documentation.

Reference
-------------
Please cite

> Widmann, A., Schröger, E., & Maess, B. (2015). Digital filter design for electrophysiological data - a practical approach. Journal of Neuroscience Methods, 250, 34-46.

if you have used functions from the EEGLAB firfilt plugin in your manuscript.

Version History
---------------
v2.8 - Added usefftfilt option to pop_eegfiltnew()

v2.7 - handle better boundary events

v2.7.1 - better handling of boundary events
