from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0059_report_gnomad_coverage_data_type'),
    ]

    operations = [
        migrations.RunSQL(
            sql="""
                ALTER TABLE variant_exlovd
                    DROP COLUMN IF EXISTS "Exon",
                    DROP COLUMN IF EXISTS "DNA_Change",
                    DROP COLUMN IF EXISTS "BIC_DNA_Change",
                    DROP COLUMN IF EXISTS "DBID",
                    DROP COLUMN IF EXISTS "Product_Of_LRs",
                    DROP COLUMN IF EXISTS "Protein_Change",
                    DROP COLUMN IF EXISTS "Source_URL",
                    DROP COLUMN IF EXISTS "Splicing_Prior_P";
            """,
            reverse_sql="""
                ALTER TABLE variant_exlovd
                    ADD COLUMN IF NOT EXISTS "Exon"           text NOT NULL DEFAULT '-',
                    ADD COLUMN IF NOT EXISTS "DNA_Change"     text NOT NULL DEFAULT '-',
                    ADD COLUMN IF NOT EXISTS "BIC_DNA_Change" text NOT NULL DEFAULT '-',
                    ADD COLUMN IF NOT EXISTS "DBID"           text NOT NULL DEFAULT '-',
                    ADD COLUMN IF NOT EXISTS "Product_Of_LRs" text NOT NULL DEFAULT '-',
                    ADD COLUMN IF NOT EXISTS "Protein_Change"   text NOT NULL DEFAULT '-',
                    ADD COLUMN IF NOT EXISTS "Source_URL"       text NOT NULL DEFAULT '-',
                    ADD COLUMN IF NOT EXISTS "Splicing_Prior_P" text NOT NULL DEFAULT '-';
            """,
        ),
    ]
