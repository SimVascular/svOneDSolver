#!/bin/bash
for f in *.h; do
	cat header $f > $f.tmp
	mv $f.tmp $f
done
