KEYS_FILE="/tmp/zotero_probe_docx_CgIJHw/zotero_item_keys.txt"
DB_COPY="/tmp/zotero_probe_db_C69iS0.sqlite"

# If the DB copy was cleaned up, re-copy it
if [ ! -f "$DB_COPY" ]; then
    cp "/mnt/c/Users/Luised94/Zotero/zotero.sqlite" "$DB_COPY"
fi

IN_CLAUSE=$(sed "s/.*/'&'/" "$KEYS_FILE" | paste -sd,)

echo "=== Keys with NO linkMode=2 PDF attachment ==="
sqlite3 -header -column "$DB_COPY" \
"SELECT pi.key,
        COALESCE(it.typeName, 'unknown') AS itemType,
        (SELECT GROUP_CONCAT(iv.value, ' // ')
         FROM itemData id2
         JOIN itemDataValues iv ON id2.valueID = iv.valueID
         JOIN fields f ON id2.fieldID = f.fieldID
         WHERE id2.itemID = pi.itemID AND f.fieldName = 'title'
        ) AS title
 FROM items pi
 LEFT JOIN itemTypes it ON pi.itemTypeID = it.itemTypeID
 WHERE pi.key IN (${IN_CLAUSE})
   AND pi.key NOT IN (
       SELECT pi2.key
       FROM items pi2
       JOIN itemAttachments ia ON ia.parentItemID = pi2.itemID
       WHERE ia.linkMode = 2 AND ia.contentType = 'application/pdf'
   );"

echo ""
echo "=== What attachments DO these keys have (if any)? ==="
sqlite3 -header -column "$DB_COPY" \
"SELECT pi.key,
        ia.linkMode,
        ia.contentType,
        ia.path
 FROM items pi
 LEFT JOIN itemAttachments ia ON ia.parentItemID = pi.itemID
 WHERE pi.key IN (${IN_CLAUSE})
   AND pi.key NOT IN (
       SELECT pi2.key
       FROM items pi2
       JOIN itemAttachments ia2 ON ia2.parentItemID = pi2.itemID
       WHERE ia2.linkMode = 2 AND ia2.contentType = 'application/pdf'
   );"
