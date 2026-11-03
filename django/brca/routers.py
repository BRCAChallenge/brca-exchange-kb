
# Split the community database tables into a separate named schema
# from the brca data, so that it can roll over between data releases
# instead of having to be re-added with each release.
# This router controls only the apps listed in route_app_labels
# It has no opinion (returns None) for 'data' app

class AppsRouter:
	"""
	A router to control all database operations on models in
	all applications other than "data"; ie, all applications
	where info persists in the db rather than being updated
	by running the pipeline.
	"""
	route_app_labels = {
        "admin",
        "auth",
        "authtoken",
        "contenttypes",
        "sessions",
        "sites",
        "users",
	}

	apps_db_alias = "community" # Which db to use?

	def db_for_read(self, model, **hints):
		if model._meta.app_label in self.route_app_labels:
			return self.apps_db_alias
		return None

	def db_for_write(self, model, **hints):
		if model._meta.app_label in self.route_app_labels:
			return self.apps_db_alias
		return None

	def allow_migrate(self, db, app_label, model_name=None, **hints):
		"""
		Allow relevant apps to migrate only if the requested db is our chosen
		community db. No opinion for other apps.
		"""
		if app_label in self.route_app_labels:
			return db == self.apps_db_alias
		return None
