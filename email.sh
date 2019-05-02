#!/usr/bin/env bash

SUBJECT="Calculations are finished"
TO="c.blair@student.rug.nl"
MESSAGE="message.txt"

/usr/bin/mail -s "$SUBJECT" "$TO" < $MESSAGE