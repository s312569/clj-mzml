# clj-mzml

A Clojure library for parsing mzML files.

## Usage

Install as:

```clj
[clj-mzml "0.1.0"]
```

Require:

```clj
(:require [clj-mzml.core :as mz])
```

Returns a hash map of parameters and other information when
`mzml-parameters` is called on an mzML file:

```clj
clj-mzml.core> (def tf "/path/to/file.mzML")
#'clj-mzml.core/tf
clj-mzml.core> (mzml-parameters tf)
{:samples {}, :file-content {:cvparams {"MS:1000579" ...
```

Keys in the returned hash are as follows:

```clj
clj-mzml.core> (keys (mzml-parameters tf))
(:samples :file-content :instrument-config :contacts :source-files :cvlist
:spectrum-count :run-parameter :default-data-proc-ref :scan-settings
:referenceable :software :data-processing)
clj-mzml.core> 
```

To access spectra use `spectra-seq` on a buffered reader opened on an
mzML file. Returns a lazy list of spectra:

```clj
clj-mzml.core> (with-open [r (io/reader tf)]
                 (->> (spectra-seq r)
                      first))
{:index "0", :id "sample=1 period=1 cycle=1 experiment=1", ...
```

Keys in the returned hashes are as follows:

```clj
clj-mzml.core> (with-open [r (io/reader tf)]
                 (->> (spectra-seq r)
                      first
                      keys))
(:index :id :defaultArrayLength :cvparams :references
:scans :precursors :products :binaries)
clj-mzml.core>
```

## License

Copyright © 2016 Jason Mulvenna

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
