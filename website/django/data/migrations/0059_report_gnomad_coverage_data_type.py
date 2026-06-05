# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0058_add_provisional_evidence_fields_to_variant'),
    ]

    operations = [
        migrations.RunSQL(
            sql=[
                "ALTER TABLE report_gnomad ADD COLUMN IF NOT EXISTS coverage TEXT NOT NULL DEFAULT '-'",
                "ALTER TABLE report_gnomad ADD COLUMN IF NOT EXISTS data_type TEXT NOT NULL DEFAULT '-'",
                'ALTER TABLE report_gnomad DROP CONSTRAINT IF EXISTS "report_gnomad_VRS_Digest_id_version_key"',
                'ALTER TABLE report_gnomad ADD CONSTRAINT "report_gnomad_VRS_Digest_id_version_data_type_key" UNIQUE ("VRS_Digest_id", version, data_type)',
            ],
            reverse_sql=[
                'ALTER TABLE report_gnomad DROP CONSTRAINT IF EXISTS "report_gnomad_VRS_Digest_id_version_data_type_key"',
                'ALTER TABLE report_gnomad ADD CONSTRAINT "report_gnomad_VRS_Digest_id_version_key" UNIQUE ("VRS_Digest_id", version)',
                "ALTER TABLE report_gnomad DROP COLUMN IF EXISTS data_type",
                "ALTER TABLE report_gnomad DROP COLUMN IF EXISTS coverage",
            ],
        ),
        # Keep Django's state in sync with the database schema.
        migrations.AddField(
            model_name='report_in_gnomad',
            name='coverage',
            field=models.TextField(default='-'),
        ),
        migrations.AddField(
            model_name='report_in_gnomad',
            name='data_type',
            field=models.TextField(default='-'),
        ),
        migrations.AlterUniqueTogether(
            name='report_in_gnomad',
            unique_together={('VRS_Digest', 'version', 'data_type')},
        ),
    ]
