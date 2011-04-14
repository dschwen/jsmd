all: jsmd.min.js

jsmd.min.js: jsmd.js
	yui-compressor --type js jsmd.js > jsmd.min.js

