#!/usr/bin/env bash
# Update the Gandi LiveDNS A record for text2gene.org to the current external IP.
# Run via cron every 15 minutes to handle dynamic IP changes from ISP.
#
# Crontab entry:
#   */15 * * * * /opt/text2gene2/deploy/update-dns.sh >> /var/log/text2gene-dns.log 2>&1

set -euo pipefail

DOMAIN="text2gene.org"
GANDI_KEY="5afd868472bea63410413fa4983531b2bc473a49"
CACHE_FILE="/tmp/text2gene_last_ip"

CURRENT_IP=$(curl -s https://api.ipify.org)

if [ -f "$CACHE_FILE" ] && [ "$(cat "$CACHE_FILE")" = "$CURRENT_IP" ]; then
    exit 0   # IP hasn't changed, nothing to do
fi

echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) Updating $DOMAIN A record to $CURRENT_IP"

# Update @ (apex) record
curl -s -X PUT "https://api.gandi.net/v5/livedns/domains/$DOMAIN/records/@/A" \
  -H "Authorization: Bearer $GANDI_KEY" \
  -H "Content-Type: application/json" \
  -d "{\"rrset_values\": [\"$CURRENT_IP\"], \"rrset_ttl\": 300}"

# Update www record
curl -s -X PUT "https://api.gandi.net/v5/livedns/domains/$DOMAIN/records/www/A" \
  -H "Authorization: Bearer $GANDI_KEY" \
  -H "Content-Type: application/json" \
  -d "{\"rrset_values\": [\"$CURRENT_IP\"], \"rrset_ttl\": 300}"

echo "$CURRENT_IP" > "$CACHE_FILE"
echo "Done."
